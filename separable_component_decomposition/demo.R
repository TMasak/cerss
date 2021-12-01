# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("library_PIP.R")

Ns <- 2^(7:15)
ERR <- array(0,c(length(Ns),4)) # for storing results

K1 <- 50
K2 <- 50
beta <- 0.7

# for the purposes of data generation and error calculation - runs a minute or two
C <- get_gneiting_cov(K1,K2,beta) 
C_mat <- tensor2matrix(C)
EIG <- eigen(C_mat)
C_half <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)

nn=2
# for(nn in 1:length(Ns)){
  N <- Ns[nn]
  # generate data
  X <- array(0,c(N,K1,K2))
  set.seed(517)
  for(n in 1:N) X[n,,] <- C_half %*% rnorm(K1*K2)
  # estimate first 3 terms in separable expansion
  Res <- separable_expansion(X,3)
  
  # explicit reconstruction of the estimators for the purposes of error calculation
  C_hat <- Res$sigma[1]*aperm(outer(Res$A[1,,],Res$B[1,,]),c(1,3,2,4))
  ERR[nn,1] <- frobenius(C-C_hat)/frobenius(C)
  C_hat <- C_hat + Res$sigma[2]*aperm(outer(Res$A[2,,],Res$B[2,,]),c(1,3,2,4))
  ERR[nn,2] <- frobenius(C-C_hat)/frobenius(C)
  C_hat <- C_hat + Res$sigma[3]*aperm(outer(Res$A[3,,],Res$B[3,,]),c(1,3,2,4))
  ERR[nn,3] <- frobenius(C-C_hat)/frobenius(C)
  
  # what R is chosen by cross-validation - runs a minute or two
  ERR[nn,4] <- CV2(X)
# }

ERR # errors for R=1,2,3 (first three columns) and which R was chosen via cross-validation (final column)
N   # at this sample size, cf. Figure 1 (left)


