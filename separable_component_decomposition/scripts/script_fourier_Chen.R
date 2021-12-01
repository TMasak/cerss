source("library_PIP.R")
library(tensorA)

estimate_separable <- function(Y, W1=NULL, W2=NULL, max_step=1, delta=0, lambda=NULL, standardize=TRUE,
                               cov_1=NULL, cov_2=NULL)
  # INPUT:     Y - array of order N x K1 x K2 representing N observations of separable process observed on K1 x K2 grid
  #           W1 - initial guess for temporal covariance (delta-shifted identity by default)
  #           W2 - initial guess for spatial covariance (delta-shifted identity by default)
  #     max_step - maximum number of iterations
  #        delta - discrete size of the corrupted band (delta=1 corresponds to corrupted diagonal)
  #       lambda - if provided, inverse weights regularized by lambda are used (lambda=0 corresponds to MLE)
  #  standardize - whether we should scale the result. FALSE makes sense only for subsequent matrix completion
  # cov_1, cov_2 - true covariances, only for for calculating errors throughout the iterations
{
  N <- dim(Y)[1]
  K1 <- dim(Y)[2]
  K2 <- dim(Y)[3]
  if(length(W1)==0){ W1 <- matrix(0,K1,K1) 
  indx <- 1:(K1-delta)
  W1[cbind(indx+delta,indx)] <- W1[cbind(indx,indx+delta)] <- 1 }
  if(length(W2)==0){ W2 <- matrix(0,K2,K2) 
  indx <- 1:(K2-delta)
  W2[cbind(indx+delta,indx)] <- W2[cbind(indx,indx+delta)] <- 1 }
  Ycenter <- sweep(Y, c(2,3), apply(Y, c(2,3), mean))
  
  iters <- 0
  ERROR <- rep(0,max_step)
  while(iters < max_step){
    covariance_1 <- array(0, dim = c(K1,K1))
    covariance_2 <- array(0, dim = c(K2,K2))
    iters <- iters+1
    for(n in 1:N){
      # if(n %% 10 == 0){ print(n) }
      covariance_1 <- covariance_1 + ( Ycenter[n,,] %*% W2 %*% t(Ycenter[n,,]) )
      covariance_2 <- covariance_2 + ( t(Ycenter[n,,]) %*% W1 %*% Ycenter[n,,] )
    }
    covariance_1 <- covariance_1/N
    covariance_2 <- covariance_2/N
    # symmetrize and possibly flip signs (if shifted traces negative)
    covariance_1 <- (covariance_1 + t(covariance_1))/2
    covariance_2 <- (covariance_2 + t(covariance_2))/2
    covariance_1 <- covariance_1*sign(sum(diag(covariance_1)))
    covariance_2 <- covariance_2*sign(sum(diag(covariance_2)))
    if(standardize){ 
      covariance_1 <- covariance_1/sqrt(abs(sum(covariance_1*W2)))
      covariance_2 <- covariance_2/sqrt(abs(sum(covariance_2*W1))) 
    }
    W1 <- covariance_1
    W2 <- covariance_2
    if(length(lambda)>0){
      W1 <- solve(W1 + diag(lambda,K1))
      W2 <- solve(W2 + diag(lambda,K2))
    }
    if(length(cov_1)>0) { ERROR[iters] <- frobenius_separable(cov_1,cov_2,covariance_1,covariance_2)/
      frobenius_separable(cov_1,cov_2)                            }
  }
  return(list(cov1 = covariance_1, cov2 = covariance_2, err = ERROR))
}

calculate_covariance_tensor <- function(X){
  ### Calculate covariance of a 2D process.
  ### Input:  X - 2D process in form X[n, x, y], for
  #           n=1,..,N (no. of observations), x,y=1,..,K (grid size)
  ### Output: covariance of X as 4D array
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  X <- to.tensor(X)
  covariance <- var.tensor(X,along="I1")
  covariance <- array(as.vector(covariance), c(K1,K2,K1,K2))
  return(covariance)
}

estimate_Chen <- function(X){
  C_emp <- calculate_covariance_tensor(X)
  Est <- estimate_separable(X,delta=0)
  A1 <- Est$cov1
  A2 <- Est$cov2
  EIG1 <- eigen(A1)
  EIG2 <- eigen(A2)
  I <- min(which(cumsum(EIG1$values)/sum(EIG1$values)>0.95)) #I <- 50
  J <- min(which(cumsum(EIG2$values)/sum(EIG2$values)>0.95)) #J <- 50
  gamma <- array(0,c(I,J))
  Chen <- array(0,dim(C_emp))
  for(k in 1:I){
    for(l in 1:J){
      pom <- outer(EIG1$vectors[,k],EIG2$vectors[,l])
      pom <- outer(pom,pom)
      gamma[k,l] <- sum(C_emp*pom)
      Chen <- Chen + gamma[k,l]*pom
    }
  }
  return(Chen)
}

spocitej <- function(sim){
  Ns <- 2^(4:12)
  ERR <- array(0,c(length(Ns),25)) # sample size x [R=1,2,3,4 + empirical + 
  # prediction times (for R=1,2,3,4 + empirical) + esimation time (for R=4 + empirical) + 3 empty +
  # R=1,2,3,4 + empirical with different denominator]
  rate_out <- 2
  rate_in <- 1.3
  
  K <- K1 <- K2 <- 50
  R <- 3
  
  As <- array(0,c(R,K,K))
  Bs <- array(0,c(R,K,K))
  Sigmas <- rep(0,K)
  A_halfs <- array(0,c(R,K,K))
  B_halfs <- array(0,c(R,K,K))
  
  V <- array(0,c(K,K))
  x <- (1:K)/(K+1)
  for(k in 1:K){
    if(k %% 2 == 1) V[,k] <- sin(2*pi*x*(k %/% 2 + 1)) else V[,k] <- cos(2*pi*x*(k %/% 2))
    V[,k] <- V[,k]/frobenius(V[,k])
  }
  order1 <- c(seq(1,K,by=R),seq(2,K,by=R),seq(3,K,by=R))
  order2 <- c(seq(2,K,by=R),seq(3,K,by=R),seq(1,K,by=R))
  order3 <- c(seq(3,K,by=R),seq(1,K,by=R),seq(2,K,by=R))
  As[1,,] <- V[,order1] %*% diag(rate_in^(K:1)) %*% t(V[,order1])
  As[2,,] <- V[,order2] %*% diag(rate_in^(K:1)) %*% t(V[,order2])
  As[3,,] <- V[,order3] %*% diag(rate_in^(K:1)) %*% t(V[,order3])
  As[1,,] <- As[1,,]/frobenius(As[1,,])
  As[2,,] <- As[2,,]/frobenius(As[2,,])
  As[3,,] <- As[3,,]/frobenius(As[3,,])
  As[1,,] <- As[1,,] + array(0.002,c(K,K))
  As[2,,] <- As[2,,] + array(0.002,c(K,K))
  As[3,,] <- As[3,,] + array(0.002,c(K,K))
  As[1,,] <- As[1,,]/frobenius(As[1,,])
  As[2,,] <- As[2,,]/frobenius(As[2,,])
  As[3,,] <- As[3,,]/frobenius(As[3,,])
  
  sum(As[1,,]*As[2,,])
  sum(As[1,,]*As[3,,])
  sum(As[2,,]*As[3,,])

  sigmas <- rate_out^(R:1)
  Bs[1,,] <- sigmas[1]*As[2,,]
  Bs[2,,] <- sigmas[2]*As[3,,]
  Bs[3,,] <- sigmas[3]*As[1,,]
  C <- array(0,c(K1,K2,K1,K2))
  for(r in 1:R){
    C <- C + aperm(outer(As[r,,],Bs[r,,]),c(1,3,2,4))
    EIG <- eigen(As[r,,])
    A_halfs[r,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
    EIG <- eigen(Bs[r,,])
    B_halfs[r,,] <- EIG$vectors %*% diag(sqrt(pmax(EIG$values,0))) %*% t(EIG$vectors)
  }
  
  for(nn in 1:length(Ns)){
    print(c("Sim",sim,"nn",nn))
    # print(Sys.time())
    N <- Ns[nn]
    akt_err <- array(0,c(10,100)) # [R=1,2,3,4 + empirical] x [test samples]
    akt_times <- array(0,c(5,100))
    set.seed(517*sim)
    X <- array(0,c(N,K,K))
    Xtest <- array(0,c(100,K1,K2))
    for(n in 1:N){
      for(r in 1:R){
        X[n,,] <- X[n,,] + A_halfs[r,,] %*% matrix(rnorm(K*K),ncol=K) %*% B_halfs[r,,]
      }
    }
    for(n in 1:100){
      for(r in 1:R){
        Xtest[n,,] <- Xtest[n,,] + A_halfs[r,,] %*% matrix(rnorm(K*K),ncol=K) %*% B_halfs[r,,]
      }
    }
    # X <- X + 1e-2*rnorm(N*K1*K2)
    # Xtest <- Xtest + 1e-2*rnorm(100*K1*K2)
    
    C_hat <- estimate_Chen(X)
    ERR[nn,14] <- frobenius(C-C_hat)/frobenius(C)
    C_mat <- tensor2matrix(C_hat)
    C_mat <- C_mat + 1e-2*diag(K1*K2)

    I <- J <- 50
    Ind <- array(1,c(50,50))
    Ind[I,] <- 0
    Ind[,J] <- 0

    # predition for empirical missing for now
    C22 <- C_mat[c(Ind)==1,c(Ind)==1]
    C12 <- C_mat[c(Ind)==0,c(Ind)==1]
    for(n in 1:100){
      t <- Sys.time()
      Xhat <- Xtest[n,,]
      Xhat[Ind==0] <- C12 %*% solve(C22,c(Xtest[n,-I,-J]))
      akt_times[5,n] <- difftime(Sys.time(),t,units="secs")
      akt_err[5,n] <- frobenius(Xtest[n,,] - Xhat[,])/frobenius(Xtest[n,,])
      akt_err[10,n] <- frobenius(Xtest[n,,] - Xhat[,])/frobenius(Xtest[n,,]*(!Ind))
    }

    ERR[nn,16:25] <- rowMeans(akt_err)
    ERR[nn,6:10] <- rowMeans(akt_times)
  }
  return(ERR)
}

# Data <- spocitej(1)

library(parallel)
Data <- mclapply(1:25,spocitej,mc.cores=25)
save(Data, file="data_fourier_Chen.RData")


