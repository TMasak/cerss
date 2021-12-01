# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("library_PIP.R")

spocitej <- function(sim){
  Ns <- 2^(7:15)
  ERR <- array(0,c(length(Ns),20)) # sample size x [R=1,2,3,4 + empirical + 
  # prediction times (for R=1,2,3,4 + empirical) + esimation time (for R=4 + empirical) + 3 empty +
  # R=1,2,3,4 + empirical with different denominator]
  K1 <- 50
  K2 <- 50
  beta <- 0.7
  
  C <- get_gneiting_cov(K1,K2,beta) # for the purposes of error calculation
  C_mat <- tensor2matrix(C)
  EIG <- eigen(C_mat)
  C_half <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
  
  for(nn in 1:length(Ns)){
    print(c("Sim",sim,"nn",nn))
    akt_err <- array(0,c(10,100)) # [R=1,2,3,4 + empirical] x [test samples]
    akt_times <- array(0,c(5,100))
    N <- Ns[nn]
    X <- array(0,c(N,K1,K2))
    Xtest <- array(0,c(100,K1,K2))
    set.seed(517*sim)
    for(n in 1:N) X[n,,] <- C_half %*% rnorm(K1*K2)
    for(n in 1:100) Xtest[n,,] <- C_half %*% rnorm(K1*K2)
    # X <- X + 1e-2*rnorm(N*K1*K2)
    # Xtest <- Xtest + 5e-2*rnorm(100*K1*K2)
    
    t <- Sys.time() 
    Res <- separable_expansion(X,3,10)
    ERR[nn,11] <- difftime(Sys.time(),t,units="secs")
    
    t <- Sys.time()
    C_hat <- array(0, c(K1,K2,K1,K2))
    for(n in 1:N){
      if(n %% 10 ==0) print(n)
      C_hat <- C_hat + outer(X[n,,],X[n,,])
    }
    C_hat <- C_hat/N
    ERR[nn,12] <- difftime(Sys.time(),t,units="secs")
    C_mat <- tensor2matrix(C_hat)
    C_mat <- C_mat + 1e-2*diag(K1*K2)
    
    # I <- sample(1:K1,floor(K1/3))
    # J <- sample(1:K2,floor(K2/3))
    I <- J <- 50
    Ind <- array(1,c(50,50))
    Ind[I,] <- 0
    Ind[,J] <- 0
    for(r in 1:3){
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
        akt_err[r,n] <- frobenius(Xtest[n,,] - Xhat)/frobenius(Xtest[n,,])
        akt_err[r+5,n] <- frobenius(Xtest[n,,] - Xhat[,])/frobenius(Xtest[n,,]*(!Ind))
      }
    }
    
    
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
    
    ERR[nn,c(1:5,16:20)] <- rowMeans(akt_err)
    ERR[nn,6:10] <- rowMeans(akt_times)
  }
  return(ERR)
}

library(parallel)
Data <- mclapply(1:25,spocitej,mc.cores=25)
save(Data, file="data_gneiting_prediction.RData")


