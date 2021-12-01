source("library_PIP.R")

spocitej <- function(sim){
  Ns <- 2^(7:15)
  ERR <- array(0,c(length(Ns),16,20)) # parameters x [7 methods (R=1,..,7) + 7 singular values + 2x CV] x 20 simulation runs
  
  K1 <- 50
  K2 <- 50
  beta <- 0.7
  
  C <- get_gneiting_cov(K1,K2,beta) # for the purposes of error calculation
  C_mat <- tensor2matrix(C)
  EIG <- eigen(C_mat)
  C_half <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
  
  for(nn in 1:length(Ns)){
    print(c("Sim",sim,"nn",nn))
    print(Sys.time())
    N <- Ns[nn]
    X <- array(0,c(N,K1,K2))
    set.seed(517*sim)
    for(n in 1:N) X[n,,] <- C_half %*% rnorm(K1*K2)
    Res <- separable_expansion(X,7,10)
    
    C_hat <- Res$sigma[1]*aperm(outer(Res$A[1,,],Res$B[1,,]),c(1,3,2,4))
    ERR[nn,1,sim] <- frobenius(C-C_hat)/frobenius(C)
    C_hat <- C_hat + Res$sigma[2]*aperm(outer(Res$A[2,,],Res$B[2,,]),c(1,3,2,4))
    ERR[nn,2,sim] <- frobenius(C-C_hat)/frobenius(C)
    C_hat <- C_hat + Res$sigma[3]*aperm(outer(Res$A[3,,],Res$B[3,,]),c(1,3,2,4))
    ERR[nn,3,sim] <- frobenius(C-C_hat)/frobenius(C)
    C_hat <- C_hat + Res$sigma[4]*aperm(outer(Res$A[4,,],Res$B[4,,]),c(1,3,2,4))
    ERR[nn,4,sim] <- frobenius(C-C_hat)/frobenius(C)
    
    C_hat <- C_hat + Res$sigma[5]*aperm(outer(Res$A[5,,],Res$B[5,,]),c(1,3,2,4))
    ERR[nn,5,sim] <- frobenius(C-C_hat)/frobenius(C)
    C_hat <- C_hat + Res$sigma[6]*aperm(outer(Res$A[6,,],Res$B[6,,]),c(1,3,2,4))
    ERR[nn,6,sim] <- frobenius(C-C_hat)/frobenius(C)
    C_hat <- C_hat + Res$sigma[7]*aperm(outer(Res$A[7,,],Res$B[7,,]),c(1,3,2,4))
    ERR[nn,7,sim] <- frobenius(C-C_hat)/frobenius(C)
    
    ERR[nn,8:14,sim] <- Res$sigma
    # ERR[nn,15,sim] <- CV(X)
    ERR[nn,16,sim] <- CV2(X)
  }
  return(ERR)
}

library(parallel)
Data <- mclapply(1:25,spocitej,mc.cores=25)
save(Data, file="data_gneiting_estimation.RData")


