# setwd("C:/Users/Tomas/Documents/Skola/EPFL/R/sparse_separable/2020_11_27")
source("library_smooth.R")
source("library_kronPCA.R")

simulate <- function(sim){
  N <- 100       # no. of curves
  K1 <- 20; K2 <- 20 # grid sizes 
  ERR <- array(0,c(7,7)) # iter=1,2,3, runtime(for 2), bandwidths , fully observed x different percentages observed
  Perc <- c(1,2,5,10,20,40,70)
  A <- get_cov(K1,"brownian")
  A <- A/frobenius(A)
  A <- A+array(0.2/K1^2,c(K1,K1))
  A <- A/frobenius(A)
  B <- get_cov(K2,"brownian")
  B <- B/frobenius(B)
  B <- B+array(0.2/K2^2,c(K2,K2))
  B <- B/frobenius(B)
  for(pp in 1:7){ 
    perc <- Perc[pp]
    print(perc)
    set.seed(17*sim)
    X <- simulate_data(N,perc,A,B,WN_flag=F)
    X <- X + 0.05*array(rnorm(N*K1*K2),c(N,K1,K2))
    
    Res <- estimate_sparse_sim(X,WN_flag=T,maxiter=3,Cov1=A,Cov2=B)
    ERR[1:3,pp] <- Res$err
    
    ERR[5,pp] <- Res$bw1
    ERR[6,pp] <- Res$bw2
    
    if(pp < 5){
      C_hat <- smooth4D(X,Res$bw1,Res$bw2,F)
      ERR[4,pp] <- frobenius( C_hat$C - aperm(outer(A,B),c(1,3,2,4)) )
    }
    
  }
  X <- simulate_data(N,100,A,B,WN_flag=F)
  Res <- kronPCA(X,1,10)
  ERR[7,7] <- frobenius( Res$sigma[1]*outer(Res$A[1,,],Res$B[1,,]) - outer(A,B) )
  
  return(ERR)
}

library(parallel)
Data <- mclapply(1:100,simulate,mc.cores=25)
save(Data, file="sparse_brownian_data.RData")
