# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("library.R")
library(RSpectra)

VanLoansPerm <- function(C_mat,K1,K2){
  return( tensor2matrix(aperm(matrix2tensor(C_mat,K1,K2),c(1,3,2,4))) )
}

tensor2matrix <- function(C){
  # transforms a covariance tensor into the proper covariance matrix
  K1 <- dim(C)[1]
  K2 <- dim(C)[2]
  C_mat <- matrix(c(C),ncol=K1*K2)
  return(C_mat)
}

matrix2tensor <- function(C_mat,K1,K2){
  # transforms a covariance matrix into the proper covariance tensor, dimensions must be provided
  C <- array(c(C_mat),c(K1,K2,K1,K2))
  return(C)
}

sim <- 1
RES <- array(0,c(11,12))
K <- 100
for(buf in 0:9){
  print(buf)
  delta <- 2*buf+1
  Dat <- create_data(seed=17*sim,buf=buf,A1="brownian",A2="brownian",mask="ar",tau=0.5,K=K,snr=3)
  
  B <- expand_b(Dat$B,delta,K)
  C <- aperm(outer(Dat$A1, Dat$A2),c(1,3,2,4)) + B
  C <- VanLoansPerm(tensor2matrix(C),K,K)
  SVD <- svds(C,1,0,0)
  RES[buf+2,1] <- SVD$d[1]/frobenius(C)
}
for(buf in 0:9){
  print(buf)
  delta <- 2*buf+1
  Dat <- create_data(seed=17*sim,buf=buf,A1="brownian",A2="brownian",mask="triangular",tau=0.5,K=K,snr=3)
  
  B <- expand_b(Dat$B,delta,K)
  C <- aperm(outer(Dat$A1, Dat$A2),c(1,3,2,4)) + B
  C <- VanLoansPerm(tensor2matrix(C),K,K)
  SVD <- svds(C,1,0,0)
  RES[buf+2,2] <- SVD$d[1]/frobenius(C)
}
for(buf in 0:9){
  print(buf)
  delta <- 2*buf+1
  Dat <- create_data(seed=17*sim,buf=buf,A1="legendre",A2="legendre",mask="ar",tau=0.5,K=K,snr=3)
  
  B <- expand_b(Dat$B,delta,K)
  C <- aperm(outer(Dat$A1, Dat$A2),c(1,3,2,4)) + B
  C <- VanLoansPerm(tensor2matrix(C),K,K)
  SVD <- svds(C,1,0,0)
  RES[buf+2,3] <- SVD$d[1]/frobenius(C)
}
for(buf in 0:9){
  print(buf)
  delta <- 2*buf+1
  Dat <- create_data(seed=17*sim,buf=buf,A1="legendre",A2="legendre",mask="triangular",tau=0.5,K=K,snr=3)
  
  B <- expand_b(Dat$B,delta,K)
  C <- aperm(outer(Dat$A1, Dat$A2),c(1,3,2,4)) + B
  C <- VanLoansPerm(tensor2matrix(C),K,K)
  SVD <- svds(C,1,0,0)
  RES[buf+2,4] <- SVD$d[1]/frobenius(C)
}
#########################################################################################################
TAU <- c(1,2,4,8,16,32,64,128,256)
buf <- 4
for(tau in 1:9){
  print(buf)
  delta <- 2*buf+1
  Dat <- create_data(seed=17*sim,buf=buf,A1="brownian",A2="brownian",mask="ar",tau=0.5,K=K,snr=TAU[tau])
  
  B <- expand_b(Dat$B,delta,K)
  C <- aperm(outer(Dat$A1, Dat$A2),c(1,3,2,4)) + B
  C <- VanLoansPerm(tensor2matrix(C),K,K)
  SVD <- svds(C,1,0,0)
  RES[tau,5] <- SVD$d[1]/frobenius(C)
}
for(tau in 1:9){
  print(buf)
  delta <- 2*buf+1
  Dat <- create_data(seed=17*sim,buf=buf,A1="brownian",A2="brownian",mask="triangular",tau=0.5,K=K,snr=TAU[tau])
  
  B <- expand_b(Dat$B,delta,K)
  C <- aperm(outer(Dat$A1, Dat$A2),c(1,3,2,4)) + B
  C <- VanLoansPerm(tensor2matrix(C),K,K)
  SVD <- svds(C,1,0,0)
  RES[tau,6] <- SVD$d[1]/frobenius(C)
}
for(tau in 1:9){
  print(buf)
  delta <- 2*buf+1
  Dat <- create_data(seed=17*sim,buf=buf,A1="legendre",A2="legendre",mask="ar",tau=0.5,K=K,snr=TAU[tau])
  
  B <- expand_b(Dat$B,delta,K)
  C <- aperm(outer(Dat$A1, Dat$A2),c(1,3,2,4)) + B
  C <- VanLoansPerm(tensor2matrix(C),K,K)
  SVD <- svds(C,1,0,0)
  RES[tau,7] <- SVD$d[1]/frobenius(C)
}
for(tau in 1:9){
  print(buf)
  delta <- 2*buf+1
  Dat <- create_data(seed=17*sim,buf=buf,A1="legendre",A2="legendre",mask="triangular",tau=0.5,K=K,snr=TAU[tau])
  
  B <- expand_b(Dat$B,delta,K)
  C <- aperm(outer(Dat$A1, Dat$A2),c(1,3,2,4)) + B
  C <- VanLoansPerm(tensor2matrix(C),K,K)
  SVD <- svds(C,1,0,0)
  RES[tau,8] <- SVD$d[1]/frobenius(C)
}
########################################################################################################
Ns <- 2^(5:11)
for(tau in 1:7){
  print(buf)
  delta <- 2*buf+1
  Dat <- create_data(seed=17*sim,buf=buf,A1="brownian",A2="brownian",mask="ar",tau=0.5,K=K,N=Ns[tau],snr=3)
  
  B <- expand_b(Dat$B,delta,K)
  C <- aperm(outer(Dat$A1, Dat$A2),c(1,3,2,4)) + B
  C <- VanLoansPerm(tensor2matrix(C),K,K)
  SVD <- svds(C,1,0,0)
  RES[tau,9] <- SVD$d[1]/frobenius(C)
}
for(tau in 1:7){
  print(buf)
  delta <- 2*buf+1
  Dat <- create_data(seed=17*sim,buf=buf,A1="brownian",A2="brownian",mask="triangular",tau=0.5,K=K,N=Ns[tau],snr=3)
  
  B <- expand_b(Dat$B,delta,K)
  C <- aperm(outer(Dat$A1, Dat$A2),c(1,3,2,4)) + B
  C <- VanLoansPerm(tensor2matrix(C),K,K)
  SVD <- svds(C,1,0,0)
  RES[tau,10] <- SVD$d[1]/frobenius(C)
}
for(tau in 1:7){
  print(buf)
  delta <- 2*buf+1
  Dat <- create_data(seed=17*sim,buf=buf,A1="legendre",A2="legendre",mask="ar",tau=0.5,K=K,N=Ns[tau],snr=3)
  
  B <- expand_b(Dat$B,delta,K)
  C <- aperm(outer(Dat$A1, Dat$A2),c(1,3,2,4)) + B
  C <- VanLoansPerm(tensor2matrix(C),K,K)
  SVD <- svds(C,1,0,0)
  RES[tau,11] <- SVD$d[1]/frobenius(C)
}
for(tau in 1:7){
  print(buf)
  delta <- 2*buf+1
  Dat <- create_data(seed=17*sim,buf=buf,A1="legendre",A2="legendre",mask="triangular",tau=0.5,K=K,N=Ns[tau],snr=3)
  
  B <- expand_b(Dat$B,delta,K)
  C <- aperm(outer(Dat$A1, Dat$A2),c(1,3,2,4)) + B
  C <- VanLoansPerm(tensor2matrix(C),K,K)
  SVD <- svds(C,1,0,0)
  RES[tau,12] <- SVD$d[1]/frobenius(C)
}
save(Norms,file="Norms.RData")
  

