source("library_PIP.R")

simulate <- function(R,eps){
  N <- 500
  K1 <- 50
  K2 <- 50
  
  As <- array(0,c(R,K1,K1))
  Bs <- array(0,c(R,K2,K2))
  A_halfs <- array(0,c(R,K1,K1))
  B_halfs <- array(0,c(R,K2,K2))
  
  As[1,,] <- brownian_cov(K1)
  As[1,,] <- As[1,,]/frobenius(As[1,,])
  Bs[1,,] <- brownian_cov(K2)
  Bs[1,,] <- Bs[1,,]/frobenius(Bs[1,,])
  EIG <- eigen(As[1,,])
  A_halfs[1,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
  EIG <- eigen(Bs[1,,])
  B_halfs[1,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
  
  for(r in 2:R){
    As[r,,] <- rnd_cov(K1)
    Bs[r,,] <- rnd_cov(K2)
    EIG <- eigen(As[r,,])
    A_halfs[r,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
    EIG <- eigen(Bs[r,,])
    B_halfs[r,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
  }
  X <- array(0,c(N,K1,K2))
  for(n in 1:N){
    for(r in 1:R){
      X[n,,] <- X[n,,] + A_halfs[r,,] %*% matrix(rnorm(K1*K2),ncol=K2) %*% B_halfs[r,,]
    }
  }
  Res <- separable_expansion(X,R,10)
  
  if(R==3) Sigma1 <- seq(0.35,0.95,by=0.1)
  if(R==5) Sigma1 <- seq(0.25,0.95,by=0.1)
  if(R==7) Sigma1 <- seq(0.15,0.95,by=0.1)
  OUT <- array(0,c(9,2))
  for(ss in 1:length(Sigma1)){
    sigma1 <- Sigma1[ss]
    Sigmas <- rep(0,R)
    Sigmas[1] <- Sigma1[ss]
    total <- Sigmas[1]
    for(r in 2:(R-1)){
      Sigmas[r] <- runif(1, max(0,1-total-(R-r)*sigma1), min(1-total,sigma1))
      total <- total + Sigmas[r]
    }
    Sigmas[R] <- 1-total
    Res$sigma <- Sigmas
    
    # target_cond_no <- 500
    # lambda <- eigvals(Res)
    # eps <- (lambda$max-target_cond_no*lambda$min)/(target_cond_no - 1)
    Res_new <- list(sigma=rep(0,R+1),A=array(0,c(R+1,K1,K1)), B=array(0,c(R+1,K2,K2)))
    Res_new$sigma[1:R] <- Res$sigma
    Res_new$A[1:R,,] <- Res$A
    Res_new$B[1:R,,] <- Res$B
    Res_new$sigma[R+1] <- eps
    Res_new$A[R+1,,] <- diag(K1)
    Res_new$B[R+1,,] <- diag(K2)
    
    Z <- matrix(rnorm(K1*K2),ncol=K2)
    Y <- apply_lhs(Res_new,Z)
    Z_hat <- pcg2(Res,Y,eps)
    
    OUT[ss,] <- c(sigma1,Z_hat$iter)
  }
  return(OUT)
}

nsim=10
Rs <- c(7,5,3)
EPSs <- c(0.5,0.2,0.1,0.05,0.02,0.01)
Data <- array(0,c(nsim,length(Rs),length(EPSs),9,2))

for(sim in 1:nsim){
  for(i in 1:length(Rs)){
    for(j in 1:length(EPSs)){
      print(c(sim,i,j))
      Data[sim,i,j,,] <- simulate(Rs[i],EPSs[j])
    }
  }
}

save(Data, file="data_inverse.RData")


