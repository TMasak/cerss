# setwd("C:/Users/Tomas/Documents/Skola/EPFL/R/shifted_partial_tracing/simulace")
source('library.R')

spocitej <- function(sim){
  RES <- array(0,c(10,9))
  for(buf in 1:6){
    print(buf)
    delta <- 2*buf+1
    K <- delta*10
    Dat <- create_data(seed=17*sim,buf=buf,A1="legendre",A2="legendre",mask="ar",tau=0.5,K=K,snr=10)
    C <- list(A1=Dat$A1, A2=Dat$A2, B=Dat$B)
    # true delta
    Chat <- estimate_all(Dat$Y,delta)
    RES[buf,1] <- frobenius_all(C,Chat)
    # partial tracing (delta=0)
    Chat <- estimate_all(Dat$Y,0)
    RES[buf,2] <- frobenius_all(C,Chat)
    # CV
    FitCV <- CV(Dat$Y,10,floor(K/4),1)
    delta_hat <- min(which.max(FitCV[1,]/FitCV[2,]))
    Chat <- estimate_all(Dat$Y,delta_hat)
    RES[buf,3] <- frobenius_all(C,Chat)
    # best sepparable approximation
    Est <- separable_expansion(Dat$Y,1)
    Chat <- list(A1=drop(Est$sigma*Est$A), A2=drop(Est$A), B=array(0,c(max(1,delta),max(1,delta))))
    RES[buf,4] <- frobenius_all(C,Chat)
    # empirical covariance
    RES[buf,5] <- frobenius_empirical(Dat$Y,Dat$A1,Dat$A2,Dat$B)
    # oracle
    # Est <- estimate_separable(Dat$X)
    # Bhat <- as.matrix(toeplitz_average(Dat$W,delta))
    # Chat <- list(A1=Est$A1, A2=Est$A2, B=Bhat)
    # RES[buf,6] <- frobenius_all(C,Chat)
    ### Inverse problem - making sure that everything is PSD
    theta <- 1e-5
    EIG <- eigen(C$A1)
    A1 <- EIG$vectors %*% diag(abs(EIG$values)) %*% t(EIG$vectors)
    EIG <- eigen(C$A2)
    A2 <- EIG$vectors %*% diag(abs(EIG$values)) %*% t(EIG$vectors)
    tev <- Re(fft(to_book_format(C$B,K)))
    tev[tev < 0] <- 0
    B <- Re(fft(tev,inverse = T)/(4*K^2))
    B <- B[1:K,1:K]
    x <- runif(K^2)
    X <- matrix(x,ncol=K)
    y <- c(A1 %*% X %*% A2) + tx(tev,x) + theta*x
    x_adi <- ADI(A1, A2, B, y, theta, rho=1e-3, adapt=T,maxiter=100,tol=10^-6)

    RES[buf,7] <- sum((x-x_adi$x)^2)/sum(x^2)
    RES[buf,8] <- x_adi$iter
    RES[buf,9] <- mean(x_adi$PCGiter)
  }
  return(RES)
}
  

library(parallel)
data <- mclapply(1:25,spocitej,mc.cores=25)
save(data, file="data_adi.RData")
