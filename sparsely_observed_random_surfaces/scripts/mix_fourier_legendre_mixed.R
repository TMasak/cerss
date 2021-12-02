# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("library_smooth.R")
source("library_kronPCA.R")

simulate <- function(sim){
  N <- 100       # no. of curves
  K1 <- 20; K2 <- 20 # grid sizes 
  ERR <- array(0,c(7,7)) # iter=1,2,3, runtime(for 2) , fully observed x different percentages observed
  Perc <- c(1,2,5,10,20,40,70)
  A <- get_cov(K1,"fourier")
  A <- A/frobenius(A)
  A <- A+array(0.2/K1^2,c(K1,K1))
  A <- A/frobenius(A)
  B <- get_cov(K2,"fourier")
  B <- B/frobenius(B)
  B <- B+array(0.2/K2^2,c(K2,K2))
  B <- B/frobenius(B)
  C1 <- aperm(outer(A,B),c(1,3,2,4))
  A <- get_cov(K1,"legendre")
  A <- A/frobenius(A)
  B <- get_cov(K2,"legendre")
  B <- B/frobenius(B)
  C2 <- aperm(outer(A,B),c(1,3,2,4))
  
  C <- 0.5*C1+0.5*C2
  C <- C/frobenius(C)
  C_mat <- tensor2matrix(C)
  EIG <- eigen(C_mat)
  C_half <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
  for(pp in 1:7){ 
    perc <- Perc[pp]
    print(perc)
    set.seed(17*sim)
    # data simulation
    X <- array(0,c(N,K1,K2))
    n_points <- floor(K1*K2*perc/100)
    for(n in 1:10){
      X_akt <- C_half %*% rnorm(K1*K2)
      inds <- sample(1:(K1*K2),K1*K2/2)
      X_akt[-inds] <- NA
      X[n,,] <- X_akt
    }
    for(n in 11:N){
      X_akt <- C_half %*% rnorm(K1*K2)
      inds <- sample(1:(K1*K2),n_points)
      X_akt[-inds] <- NA
      X[n,,] <- X_akt
    }
    X <- X + 0.05*array(rnorm(N*K1*K2),c(N,K1,K2))
    
    Res1 <- estimate_sparse_sim(X,WN_flag=T,maxiter=5,Cov1=A,Cov2=B,wopt_flag=F)
    # ERR[1:3,pp] <- Res$err
    cv1 <- CV_s(X,wopt_flag=F)
    source("library_smooth.R")
    Res2 <- estimate_sparse_sim(X,WN_flag=T,maxiter=5,Cov1=A,Cov2=B,wopt_flag=T)
    # ERR[1:3,pp] <- Res$err
    cv2 <- CV_s(X,wopt_flag=T)
    source("library_smooth2.R")
    Res3 <- estimate_sparse_sim(X,WN_flag=T,maxiter=5,Cov1=A,Cov2=B,wopt_flag=T)
    # ERR[1:3,pp] <- Res$err
    cv3 <- CV_s(X,wopt_flag=T)
    
    ERR[5,pp] <- Res$bw1
    ERR[6,pp] <- Res$bw2
    
    if(pp < 5){
      C_hat <- smooth4D(X,Res$bw1,Res$bw2,F)
      ERR[4,pp] <- frobenius( C_hat$C - C )
    }
    
  }
  set.seed(17*sim)
  X <- array(0,c(N,K1,K2))
  for(n in 1:N) X[n,,] <- C_half %*% rnorm(K1*K2)
  # X <- X + 0.05*array(rnorm(N*K1*K2),c(N,K1,K2))
  Res <- kronPCA(X,1,10)
  ERR[7,7] <- frobenius( Res$sigma[1]*aperm(outer(Res$A[1,,],Res$B[1,,]),c(1,3,2,4)) - C )
  
  return(ERR)
}

library(parallel)
Data <- mclapply(1:100,simulate,mc.cores=25)
save(Data, file="mix_fourier_legendre_data.RData")
