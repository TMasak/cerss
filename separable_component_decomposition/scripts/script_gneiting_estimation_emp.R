source("library_PIP.R")

spocitej <- function(sim){
  Ns <- 2^(7:15)
  ERR <- array(0,c(length(Ns),3,20)) # parameters x [7 methods (R=1,..,7) + 7 singular values + 2x CV] x 20 simulation runs
  
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
    Res <- separable_expansion(X,1,10)
    
    C_hat <- Res$sigma[1]*aperm(outer(Res$A[1,,],Res$B[1,,]),c(1,3,2,4))
    ERR[nn,1,sim] <- frobenius(C-C_hat)/frobenius(C)
    
    C_hat <- array(0, c(K1,K2,K1,K2))
    for(n in 1:N){
      C_hat <- C_hat + outer(X[n,,],X[n,,])
    }
    C_hat <- C_hat/N
    
    ERR[nn,2,sim] <- frobenius(C-C_hat)/frobenius(C)
    
    C_hat_mat <- tensor2matrix(C_hat)
    EIG <- eigen(C_hat_mat)
    
    r <- 15
    C_hat <- EIG$vectors[,1:r] %*% diag(EIG$values[1:r]) %*% t(EIG$vectors[,1:r])
    ERR[nn,3,sim] <- frobenius(C_mat-C_hat)/frobenius(C)
  }
  return(ERR)
}

library(parallel)
Data <- mclapply(1:25,spocitej,mc.cores=25)
save(Data, file="data_gneiting_emp.RData")


