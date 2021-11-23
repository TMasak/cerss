# setwd("C:/Users/Tomas/Documents/Skola/EPFL/R/shifted_partial_tracing")
# source("library_PIP.R")
source("library.R")

simulate <- function(sim){
  RES <- array(0,c(11,7))
    delta <- 0
    # Dat <- input[[sim]]
    Dat <- create_data(seed=17*sim,buf=0,A1="brownian",A2="brownian",mask="ar",snr=3)
    C <- list(A1=Dat$A1, A2=Dat$A2, B=array(0,c(1,1)))
    # true delta
    Chat <- estimate_all(Dat$X,delta)
    RES[1,1] <- frobenius_all(C,Chat)
    # partial tracing (delta=0)
    Chat <- estimate_all(Dat$X,0)
    RES[1,2] <- frobenius_all(C,Chat)
    # fit based CV
    FitCV <- CV(Dat$X,10,10,0)
    delta_hat <- min(which.max(FitCV[1,]/FitCV[2,]))-1
    Chat <- estimate_all(Dat$X,delta_hat)
    RES[1,3] <- frobenius_all(C,Chat)
    # new CV
    delta_hat <- min(which.min(FitCV[2,]^2 - 2*FitCV[1,]))-1
    Chat <- estimate_all(Dat$X,delta_hat)
    RES[1,4] <- frobenius_all(C,Chat)
    # best sepparable approximation
    Est <- separable_expansion(Dat$X,1)
    Chat <- list(A1=drop(Est$sigma*Est$A), A2=drop(Est$A), B=array(0,c(max(1,delta),max(1,delta))))
    RES[1,5] <- frobenius_all(C,Chat)
    # oracle
    Est <- estimate_separable(Dat$X)
    Bhat <- as.matrix(toeplitz_average(Dat$W,delta))
    Chat <- list(A1=Est$A1, A2=Est$A2, B=array(0,c(1,1)))
    RES[1,6] <- frobenius_all(C,Chat)
    # empirical covariance
    print(c("Empirical", sim))
    RES[1,7] <- frobenius_empirical(Dat$X,Dat$A1,Dat$A2,array(0,c(1,1)))
  for(buf in 0:9){
    print(buf)
    delta <- 2*buf+1
    Dat <- create_data(seed=17*sim,buf=buf,A1="brownian",A2="brownian",mask="ar",tau=0.5,snr=3)
    C <- list(A1=Dat$A1, A2=Dat$A2, B=Dat$B)
    # true delta
    Chat <- estimate_all(Dat$Y,delta)
    RES[buf+2,1] <- frobenius_all(C,Chat)
    # partial tracing (delta=0)
    Chat <- estimate_all(Dat$Y,0)
    RES[buf+2,2] <- frobenius_all(C,Chat)
    # fit based CV
    FitCV <- CV(Dat$Y,10,20,0)
    delta_hat <- min(which.max(FitCV[1,]/FitCV[2,]))-1
    Chat <- estimate_all(Dat$Y,delta_hat)
    RES[buf+2,3] <- frobenius_all(C,Chat)
    # new CV
    delta_hat <- min(which.min(FitCV[2,]^2 - 2*FitCV[1,]))-1
    Chat <- estimate_all(Dat$Y,delta_hat)
    RES[buf+2,4] <- frobenius_all(C,Chat)
    # best sepparable approximation
    Est <- separable_expansion(Dat$Y,1)
    Chat <- list(A1=drop(Est$sigma*Est$A), A2=drop(Est$A), B=array(0,c(max(1,delta),max(1,delta))))
    RES[buf+2,5] <- frobenius_all(C,Chat)
    # oracle
    Est <- estimate_separable(Dat$X)
    Bhat <- as.matrix(toeplitz_average(Dat$W,delta))
    Chat <- list(A1=Est$A1, A2=Est$A2, B=Bhat)
    RES[buf+2,6] <- frobenius_all(C,Chat)
    # empirical covariance
    print(c("Empirical", buf))
    RES[buf+2,7] <- frobenius_empirical(Dat$Y,Dat$A1,Dat$A2,Dat$B)
  }
  return(RES)
}

library(parallel)
data <- mclapply(1:25,simulate,mc.cores=25)
# data <- lapply(1:25,simulate) # comment the previous line an uncomment this one if not running on a cluster
save(data, file="brownian_ar_buf.RData")
