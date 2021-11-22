setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("library.R")

### Estimation

buf <- 4
delta <- 2*buf+1
# simulate data
Dat <- create_data(seed=17,buf=buf,A1="brownian",A2="brownian",mask="triangular",snr=3)
C <- list(A1=Dat$A1, A2=Dat$A2, B=Dat$B)
# SPT with true delta provided true delta
Chat <- estimate_all(Dat$Y,delta)
frobenius_all(C,Chat)
# partial tracing (delta=0)
Chat <- estimate_all(Dat$Y,0)
frobenius_all(C,Chat)
# CV
FitCV <- CV4(Dat$Y,10,20,1)
delta_hat <- min(which.max(FitCV[1,]/FitCV[2,]))
Chat <- estimate_all(Dat$Y,delta_hat)
frobenius_all(C,Chat)
# stability selection
Stab <- stability_selection(Dat$Y,20,1)
delta_hat <- max_prominence2(Stab[3,])
Chat <- estimate_all(Dat$Y,delta_hat)
frobenius_all(C,Chat)
# empirical covariance
frobenius_empirical(Dat$Y,Dat$A1,Dat$A2,Dat$B)
# oracle
Est <- estimate_separable(Dat$X)
Bhat <- as.matrix(toeplitz_average(Dat$W,delta))
Chat <- list(A1=Est$A1, A2=Est$A2, B=Bhat)
frobenius_all(C,Chat)

### Inverse problem

buf <- 4
delta <- 2*buf+1
K <- 100 # grid size is set to 100 automatically inside create_data()
Dat <- create_data(seed=17,buf=buf,A1="legendre",A2="legendre",mask="triangular",tau=0.5,K=K)
C <- estimate_all(Dat$Y,delta)

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

sum((x-x_adi$x)^2)/sum(x^2)
x_adi$iter
mean(x_adi$PCGiter)
