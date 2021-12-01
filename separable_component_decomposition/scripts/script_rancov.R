source("library_PIP.R")

spocitej <- function(sim){
  Ns <- 2^(7:13)
  Rates <- 2:10
  ERR <- array(0,c(length(Ns),length(Rates),4))
  
  K <- 50
  R <- 4
  ranks <- c(15,10,13,12)
  
  eigfuns <- matrix(rnorm(K*K),ncol=K)
  SVD <- svd(eigfuns)
  eigfuns <- SVD$u
  
  As <- array(0,c(R,K,K))
  Bs <- array(0,c(R,K,K))
  Sigmas <- rep(0,K)
  A_halfs <- array(0,c(R,K,K))
  B_halfs <- array(0,c(R,K,K))
  
  As[1,,] <- eigfuns[,1:ranks[1]] %*% diag(ranks[1]:1) %*% t(eigfuns[,1:ranks[1]])
  As[2,,] <- eigfuns[,(1+ranks[1]):(ranks[1]+ranks[2])] %*% diag(ranks[2]:1) %*% 
    t(eigfuns[,(1+ranks[1]):(ranks[1]+ranks[2])])
  As[3,,] <- eigfuns[,(1+ranks[1]+ranks[2]):(ranks[1]+ranks[2]+ranks[3])] %*% diag(ranks[3]:1) %*% 
    t(eigfuns[,(1+ranks[1]+ranks[2]):(ranks[1]+ranks[2]+ranks[3])])
  As[4,,] <- eigfuns[,(1+sum(ranks[1:3])):50] %*% diag(ranks[4]:1) %*% t(eigfuns[,(1+sum(ranks[1:3])):50])
  As[1,,] <- As[1,,]/frobenius(As[1,,])
  As[2,,] <- As[2,,]/frobenius(As[2,,])
  As[3,,] <- As[3,,]/frobenius(As[3,,])
  As[4,,] <- As[4,,]/frobenius(As[4,,])
  
  for(j in 1:length(Rates)){
    sigmas <- Rates[j]^(R:1)
    Bs[1,,] <- sigmas[1]*As[2,,]
    Bs[2,,] <- sigmas[2]*As[3,,]
    Bs[3,,] <- sigmas[3]*As[4,,]
    Bs[4,,] <- sigmas[4]*As[1,,]
    for(r in 1:R){
      EIG <- eigen(As[r,,])
      A_halfs[r,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
      EIG <- eigen(Bs[r,,])
      B_halfs[r,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
    }
    
    for(nn in 1:length(Ns)){
      print(c("Sim",sim,"nn",nn))
      print(Sys.time())
      N <- Ns[nn]
      set.seed(517*sim)
      X <- array(0,c(N,K,K))
      for(n in 1:N){
        for(r in 1:R){
          X[n,,] <- X[n,,] + A_halfs[r,,] %*% matrix(rnorm(K*K),ncol=K) %*% B_halfs[r,,]
        }
      }
      Res <- separable_expansion(X,R,10)
      
      Res_new <- Res
      Res_new$sigma <- Res_new$sigma[1]
      Res_new$A <- Res_new$A[1,,,drop=F]
      Res_new$B <- Res_new$B[1,,,drop=F]
      ERR[nn,j,1] <- frobenius_PCA(Res_new,As,Bs) 
      
      Res_new <- Res
      Res_new$sigma <- Res_new$sigma[1:2]
      Res_new$A <- Res_new$A[1:2,,]
      Res_new$B <- Res_new$B[1:2,,]
      ERR[nn,j,2] <- frobenius_PCA(Res_new,As,Bs) 
      
      Res_new <- Res
      Res_new$sigma <- Res_new$sigma[1:3]
      Res_new$A <- Res_new$A[1:3,,]
      Res_new$B <- Res_new$B[1:3,,]
      ERR[nn,j,3] <- frobenius_PCA(Res_new,As,Bs) 
      
      ERR[nn,j,4] <- frobenius_PCA(Res,As,Bs) 
    }
  }
  return(ERR)
}

library(parallel)
Data <- mclapply(1:10,spocitej,mc.cores=10)
save(Data, file="data_rancov.RData")










