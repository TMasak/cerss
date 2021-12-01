source("library_PIP.R")

spocitej <- function(sim){
  Ns <- 2^(4:12)
  ERR <- array(0,c(length(Ns),25)) # sample size x [R=1,2,3,4 + empirical + 
  # prediction times (for R=1,2,3,4 + empirical) + esimation time (for R=4 + empirical) + 3 empty +
  # R=1,2,3,4 + empirical with different denominator]
  rate_out <- 2
  rate_in <- 1.3
  
  K <- K1 <- K2 <- 50
  R <- 3
  
  As <- array(0,c(R,K,K))
  Bs <- array(0,c(R,K,K))
  Sigmas <- rep(0,K)
  A_halfs <- array(0,c(R,K,K))
  B_halfs <- array(0,c(R,K,K))
  
  V <- array(0,c(K,K))
  x <- (1:K)/(K+1)
  for(k in 1:K){
    if(k %% 2 == 1) V[,k] <- sin(2*pi*x*(k %/% 2 + 1)) else V[,k] <- cos(2*pi*x*(k %/% 2))
    V[,k] <- V[,k]/frobenius(V[,k])
  }
  order1 <- c(seq(1,K,by=R),seq(2,K,by=R),seq(3,K,by=R))
  order2 <- c(seq(2,K,by=R),seq(3,K,by=R),seq(1,K,by=R))
  order3 <- c(seq(3,K,by=R),seq(1,K,by=R),seq(2,K,by=R))
  As[1,,] <- V[,order1] %*% diag(rate_in^(K:1)) %*% t(V[,order1])
  As[2,,] <- V[,order2] %*% diag(rate_in^(K:1)) %*% t(V[,order2])
  As[3,,] <- V[,order3] %*% diag(rate_in^(K:1)) %*% t(V[,order3])
  As[1,,] <- As[1,,]/frobenius(As[1,,])
  As[2,,] <- As[2,,]/frobenius(As[2,,])
  As[3,,] <- As[3,,]/frobenius(As[3,,])
  As[1,,] <- As[1,,] + array(0.002,c(K,K))
  As[2,,] <- As[2,,] + array(0.002,c(K,K))
  As[3,,] <- As[3,,] + array(0.002,c(K,K))
  As[1,,] <- As[1,,]/frobenius(As[1,,])
  As[2,,] <- As[2,,]/frobenius(As[2,,])
  As[3,,] <- As[3,,]/frobenius(As[3,,])
  
  sum(As[1,,]*As[2,,])
  sum(As[1,,]*As[3,,])
  sum(As[2,,]*As[3,,])

  sigmas <- rate_out^(R:1)
  Bs[1,,] <- sigmas[1]*As[2,,]
  Bs[2,,] <- sigmas[2]*As[3,,]
  Bs[3,,] <- sigmas[3]*As[1,,]
  C <- array(0,c(K1,K2,K1,K2))
  for(r in 1:R){
    C <- C + aperm(outer(As[r,,],Bs[r,,]),c(1,3,2,4))
    EIG <- eigen(As[r,,])
    A_halfs[r,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
    EIG <- eigen(Bs[r,,])
    B_halfs[r,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
  }
  
  for(nn in 1:length(Ns)){
    print(c("Sim",sim,"nn",nn))
    # print(Sys.time())
    N <- Ns[nn]
    akt_err <- array(0,c(10,100)) # [R=1,2,3,4 + empirical] x [test samples]
    akt_times <- array(0,c(5,100))
    set.seed(517*sim)
    X <- array(0,c(N,K,K))
    Xtest <- array(0,c(100,K1,K2))
    for(n in 1:N){
      for(r in 1:R){
        X[n,,] <- X[n,,] + A_halfs[r,,] %*% matrix(rnorm(K*K),ncol=K) %*% B_halfs[r,,]
      }
    }
    for(n in 1:100){
      for(r in 1:R){
        Xtest[n,,] <- Xtest[n,,] + A_halfs[r,,] %*% matrix(rnorm(K*K),ncol=K) %*% B_halfs[r,,]
      }
    }
    # X <- X + 1e-2*rnorm(N*K1*K2)
    # Xtest <- Xtest + 1e-2*rnorm(100*K1*K2)

    t <- Sys.time()
    Res <- separable_expansion(X,R,10)
    ERR[nn,11] <- difftime(Sys.time(),t,units="secs")
    ERR[nn,13] <- CV2(X,8,3)
    
    C_hat <- Res$sigma[1]*aperm(outer(Res$A[1,,],Res$B[1,,]),c(1,3,2,4))
    ERR[nn,1] <- frobenius(C-C_hat)/frobenius(C)
    C_hat <- C_hat + Res$sigma[2]*aperm(outer(Res$A[2,,],Res$B[2,,]),c(1,3,2,4))
    ERR[nn,2] <- frobenius(C-C_hat)/frobenius(C)
    C_hat <- C_hat + Res$sigma[3]*aperm(outer(Res$A[3,,],Res$B[3,,]),c(1,3,2,4))
    ERR[nn,3] <- frobenius(C-C_hat)/frobenius(C)

    t <- Sys.time()
    C_hat <- array(0, c(K1,K2,K1,K2))
    for(n in 1:N){
      if(n %% 10 == 0) print(n)
      C_hat <- C_hat + outer(X[n,,],X[n,,])
    }
    C_hat <- C_hat/N
    ERR[nn,12] <- difftime(Sys.time(),t,units="secs")
    ERR[nn,5] <- frobenius(C-C_hat)/frobenius(C)
    C_mat <- tensor2matrix(C_hat)
    C_mat <- C_mat + 1e-2*diag(K1*K2)

    I <- J <- 50
    Ind <- array(1,c(50,50))
    Ind[I,] <- 0
    Ind[,J] <- 0
    for(r in 1:R){
      akt_Res <- Res
      akt_Res$A <- Res$A[1:r,,,drop=F]
      akt_Res$B <- Res$B[1:r,,,drop=F]
      akt_Res$sigma <- Res$sigma[1:r]
      EIG <- eigvals(akt_Res)
      if(EIG$min < 0) eps <- -EIG$min + 1e-2 else eps <- 1e-2
      for(n in 1:100){
        t <- Sys.time()
        Xhat <- BLUP(Xtest[n,,],akt_Res,I,J,eps)
        akt_times[r,n] <- difftime(Sys.time(),t,units="secs")
        akt_err[r,n] <- frobenius(Xtest[n,,] - Xhat[,])/frobenius(Xtest[n,,])
        akt_err[r+5,n] <- frobenius(Xtest[n,,] - Xhat[,])/frobenius(Xtest[n,,]*(!Ind))
      }
    }

    # predition for empirical missing for now
    C22 <- C_mat[c(Ind)==1,c(Ind)==1]
    C12 <- C_mat[c(Ind)==0,c(Ind)==1]
    for(n in 1:100){
      t <- Sys.time()
      Xhat <- Xtest[n,,]
      Xhat[Ind==0] <- C12 %*% solve(C22,c(Xtest[n,-I,-J]))
      akt_times[5,n] <- difftime(Sys.time(),t,units="secs")
      akt_err[5,n] <- frobenius(Xtest[n,,] - Xhat[,])/frobenius(Xtest[n,,])
      akt_err[10,n] <- frobenius(Xtest[n,,] - Xhat[,])/frobenius(Xtest[n,,]*(!Ind))
    }

    ERR[nn,16:25] <- rowMeans(akt_err)
    ERR[nn,6:10] <- rowMeans(akt_times)
  }
  return(ERR)
}

library(parallel)
Data <- mclapply(1:25,spocitej,mc.cores=25)
save(Data, file="data_fourier.RData")


