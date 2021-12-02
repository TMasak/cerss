# setwd("C:/Users/Tomas/Documents/Skola/EPFL/R/sparse_separable/2020_11_27")
source("library_smooth.R")
source("library_kronPCA.R")

simulate <- function(sim){
  N <- 100       # no. of curves
  K1 <- 20; K2 <- 20 # grid sizes 
  ERR <- array(0,c(7,7)) # iter=1,2,3, runtime(for 2), bandwidths , fully observed x different percentages observed
  Perc <- c(1,2,5,10,20,40,70)
  A <- get_cov(K1,"fourier")
  A <- A/frobenius(A)
  A <- A+array(0.2/K1^2,c(K1,K1))
  A <- A/frobenius(A)
  B <- get_cov(K2,"fourier")
  B <- B/frobenius(B)
  B <- B+array(0.2/K2^2,c(K2,K2))
  B <- B/frobenius(B)
  for(pp in 1:6){ 
    perc <- Perc[pp]
    print(perc)
    set.seed(17*sim)
    X <- simulate_data(N,perc,A,B,WN_flag=T)
    
    # Res <- estimate_sparse_sim(X,WN_flag=T,maxiter=3,Cov1=A,Cov2=B)
    # ERR[1:3,pp] <- Res$err
    
    t <- Sys.time()
    blah <- estimate_sparse_nonpooled(X,bwidths = c(6.7175144,6.7175144))
    ERR[4,pp] <- difftime(Sys.time(),t,units="secs")
    
    # t <- Sys.time()
    # C_hat <- smooth4D(X,3.5,3.5,F)
    # ERR[4,pp] <- frobenius( C_hat$C - aperm(outer(A,B),c(1,3,2,4)) )
    # cat(c("perc =",pp, "takes", Sys.time()-t, units(Sys.time()-t)))
    
  }
  X <- simulate_data(N,100,A,B,WN_flag=T)
  Res <- kronPCA(X,1,10)
  ERR[7,7] <- frobenius( outer(Res$A[1,,],Res$B[1,,]) - outer(A,B) )
  
  return(ERR)
}

library(parallel)
Data <- mclapply(1:25,simulate,mc.cores=25)
save(Data, file="sparse_fourier_nonpooled_data.RData")
