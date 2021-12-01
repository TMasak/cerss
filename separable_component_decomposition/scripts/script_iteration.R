source("library_PIP.R")

simulate <- function(sim){
  N <- 500
  R <- 5
  Ks <- (1:7)*30
  Kappas <- c(10,100,1000)
  OUT <- array(0, c(length(Ks),length(Kappas),2))
   
  for(i in 1:length(Ks)){
    K1 <- K2 <- Ks[i]
    As <- array(0,c(R,K1,K1))
    Bs <- array(0,c(R,K2,K2))
    A_halfs <- array(0,c(R,K1,K1))
    B_halfs <- array(0,c(R,K2,K2))
    
    As[1,,] <- brownian_cov(K1)
    As[1,,] <- As[1,,]/frobenius(As[1,,])
    Bs[1,,] <- brownian_cov(K2)
    Bs[1,,] <- Bs[1,,]/frobenius(Bs[1,,])
    EIG <- eigen(As[1,,])
    A_halfs[1,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
    EIG <- eigen(Bs[1,,])
    B_halfs[1,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
    
    set.seed(517*sim)
    for(r in 2:R){
      As[r,,] <- rnd_cov(K1)
      Bs[r,,] <- rnd_cov(K2)
      EIG <- eigen(As[r,,])
      A_halfs[r,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
      EIG <- eigen(Bs[r,,])
      B_halfs[r,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
    }
    X <- array(0,c(N,K1,K2))
    for(n in 1:N){
      for(r in 1:R){
        X[n,,] <- X[n,,] + A_halfs[r,,] %*% matrix(rnorm(K1*K2),ncol=K2) %*% B_halfs[r,,]
      }
    }
    Res <- separable_expansion(X,R,10)
    lambda <- eigvals(Res)
    
    for(j in 1:length(Kappas)){
      target_cond_no <- Kappas[j]
      Sigmas <- rep(0,R)
      Sigmas[1] <- runif(1,1/2,1)
      total <- Sigmas[1]
      for(r in 2:(R-1)){
        Sigmas[r] <- runif(1,0,1-total)
        total <- total + Sigmas[r]
      }
      Sigmas[R] <- 1-total
      Res$sigma <- Sigmas
      
      eps <- (lambda$max-target_cond_no*lambda$min)/(target_cond_no - 1)
      Res_new <- list(sigma=rep(0,R+1),A=array(0,c(R+1,K1,K1)), B=array(0,c(R+1,K2,K2)))
      Res_new$sigma[1:R] <- Res$sigma
      Res_new$A[1:R,,] <- Res$A
      Res_new$B[1:R,,] <- Res$B
      Res_new$sigma[R+1] <- eps
      Res_new$A[R+1,,] <- diag(K1)
      Res_new$B[R+1,,] <- diag(K2)
      
      Z <- matrix(rnorm(K1*K2),ncol=K2)
      Y <- apply_lhs(Res_new,Z)
      Z_hat <- pcg2(Res,Y,eps)
      OUT[i,j,1] <- frobenius(Z-Z_hat$U)/frobenius(Z)
      OUT[i,j,2] <- Z_hat$iter
    }
  }
  return(OUT)
}

library(parallel)
Data <- mclapply(1:100,simulate,mc.cores=10)
save(Data, file="data_iteration.RData")


