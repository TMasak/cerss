# setwd("C:/Users/Tomas/Documents/Skola/EPFL/R/sparse_separable/2020_11_11")
source("library_smooth.R")
source("library_kronPCA.R")

simulate <- function(sim){
  N <- 100       # no. of curves
  K1 <- 20; K2 <- 20 # grid sizes 
  ERR <- array(0,c(7,7)) # iter=1,2,3, runtime(for 2) , fully observed x different percentages observed
  Perc <- c(1,2,5,10,20,40,70)
  C <- get_gneiting_cov(K1,K2,0.7) # for the purposes of error calculation
  C <- C/frobenius(C)
  C_mat <- tensor2matrix(C)
  EIG <- eigen(C_mat)
  C_half <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
  for(pp in 2:2){ 
    perc <- Perc[pp]
    print(perc)
    set.seed(17*sim)
    # data simulation
    X <- array(0,c(N,K1,K2))
    n_points <- floor(K1*K2*perc/100)
    for(n in 1:N){
      X_akt <- C_half %*% rnorm(K1*K2)
      inds <- sample(1:(K1*K2),n_points)
      X_akt[-inds] <- NA
      X[n,,] <- X_akt
    }
    X <- X + 0.15*array(rnorm(N*K1*K2),c(N,K1,K2))
    
    # Res <- estimate_sparse_gneiting(X,WN_flag=T,maxiter=3,C=C)
    # ERR[1:3,pp] <- Res$err
    # 
    # t <- Sys.time()
    # blah <- estimate_sparse(X,bwidths = c(Res$bw1,Res$bw2))
    # ERR[4,pp] <- difftime(Sys.time(),t,units="secs")
    # 
    # ERR[5,pp] <- Res$bw1
    # ERR[6,pp] <- Res$bw2
    
    # results from above
    t <- Sys.time()
    C_hat <- smooth4D(X,NULL,NULL,T) # bandwidths manually read from commented
    ERR[7,pp] <- difftime(Sys.time(),t,units="secs")
    ERR[4,pp] <- frobenius( C_hat$C - C )
    ERR[5,pp] <- C_hat$Times[1]
    ERR[6,pp] <- C_hat$Times[2]
    
  }
  # X <- array(0,c(N,K1,K2))
  # for(n in 1:N) X[n,,] <- C_half %*% rnorm(K1*K2)
  # # X <- X + 0.05*array(rnorm(N*K1*K2),c(N,K1,K2))
  # Res <- kronPCA(X,1,10)
  # ERR[7,7] <- frobenius( Res$sigma[1]*outer(Res$A[1,,],Res$B[1,,]) - C )
  
  return(ERR)
}

library(parallel)
Data <- mclapply(1:25,simulate,mc.cores=25)
save(Data, file="sparse_gneiting_empirical_cv_data.RData")
