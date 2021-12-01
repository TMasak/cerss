require(lattice)
require(splines)

separable_expansion <- function(X,R,maxiter=10,B=NULL){
# Input: X - data array of dimensions N x K1 x K2
#        R -'kronecker' rank
#  maxiter - maximum number of iterations (to be replaced by tolerance)
#        B - array of dims R x K2 x K2, starting values for
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  if(length(B) == 0){
    B <- array(0, c(R,K2,K2))
    for(r in 1:R){
      B[r,,] <- diag(K2)
    }
  }
  A <- array(0,c(R,K1,K1))
  sigma <- rep(0,R)
  X <- sweep(X, c(2,3), apply(X, c(2,3), mean))
  for(r in 1:R){
    iter <- 0
    while(iter < maxiter){
      A[r,,] <- T1(X,A,B,sigma,r,N,K1,K2)
      A[r,,] <- A[r,,]/frobenius(A[r,,])
      B[r,,] <- T2(X,A,B,sigma,r,N,K1,K2)
      sigma[r] <- frobenius(B[r,,])
      B[r,,] <- B[r,,]/sigma[r]
      iter <- iter +1 
# if(iter>1) print(frobenius(A[r,,]-A_old)/frobenius(A[r,,]))
# A_old <- A[r,,]
    }
  }
  return(list(A=A,B=B,sigma=sigma))
}

T1 <- function (X,A,B,sigma,r,N,K1,K2){
# partial inner product w.r.t. to the first argument, i.e. returning A
  Res <- array(0,c(K1,K1))
  for(n in 1:N){
    Res <- Res + X[n,,] %*% B[r,,] %*% t(X[n,,])
  }
  Res <- Res/N
  if(r > 1){
    for(j in 1:(r-1)){
      Res <- Res - sigma[j] * sum(B[j,,]*B[r,,]) * A[j,,]
    }
  }
  return(Res)
}

T2 <- function (X,A,B,sigma,r,N,K1,K2){
# partial inner product w.r.t. to the second argument, i.e. returning B
  Res <- array(0,c(K2,K2))
  for(n in 1:N){
    Res <- Res + t(X[n,,]) %*% A[r,,] %*% X[n,,]
  }
  Res <- Res/N
  if(r>1){
    for(j in 1:(r-1)){
      Res <- Res - sigma[j] * sum(A[j,,]*A[r,,]) * B[j,,]
    }
  }
  return(Res)
}

BLUP <- function(Xnew,Res,I,J,eps=1e-2){
# Calculates the best linear unbiased predictor of the missing etries of X (currently only works for whole 
# rows with indeces in I and whole columns with indices in J) based on the covariance given in Res
  K1 <- dim(Xnew)[1]
  K2 <- dim(Xnew)[2]
  
  C22 <- Res
  C22$A <- Res$A[,-I,-I,drop=F]
  C22$B <- Res$B[,-J,-J,drop=F]
  Xobs <- Xnew[-I,-J]
  # Z <- pcg(C22,Xobs)
  Z <- pcg2(C22,Xobs,eps)
  # print(Z$iter)
  
  Zfull <- array(0,c(K1,K2))
  Zfull[-I,-J] <- Z$U
  Xhat <- apply_lhs(Res,Zfull)
  Xhat[-I,-J] <- Xobs
  
  return(list(Xhat=Xhat, iter=Z$iter))
  # C12 <- Res
  # C12$A <- Res$A[,I,-I]
  # C12$B <- Res$B[,-J,J]
  # Xpred <- apply_lhs(C12,Z$U)
  # 
}

CV <- function(X,Folds=10,maxR=7,maxiter=10){
# performs CV for the best R (between 1 and maxR) based on prediction power splitting data in a given
# number of folds
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Nnew <- (N %/% Folds)*Folds
  Ind <- matrix(sample(1:Nnew), nrow=Folds)
  ERR <- array(0,c(Folds,maxR))
  for(fold in 1:Folds){
    Xtrain <- X[-Ind[fold,],,]
    Xtest <- X[Ind[fold,],,]
    Ntest <- dim(Xtest)[1]
    
    Res <- separable_expansion(Xtrain,maxR,maxiter)
    for(r in 1:maxR){
      akt_err <- array(0,c(2,Ntest))
      akt_Res <- Res
      akt_Res$A <- Res$A[1:r,,,drop=F]
      akt_Res$B <- Res$B[1:r,,,drop=F]
      akt_Res$sigma <- Res$sigma[1:r]
      for(n in 1:Ntest){
        
        I <- sample(1:K1,floor(K1/2))
        J <- sample(1:K2,floor(K2/2))
        Xhat <- BLUP(Xtest[n,,],akt_Res,I,J)
        akt_err[1,n] <- frobenius(Xtest[n,,] - Xhat)/frobenius(Xtest[n,,])
        
        I <- sample(1:K1,floor(K1/2))
        J <- sample(1:K2,floor(K2/2))
        Xhat <- BLUP(Xtest[n,,],akt_Res,I,J)
        akt_err[2,n] <- frobenius(Xtest[n,,] - Xhat)/frobenius(Xtest[n,,])
      }
      ERR[k,r] <- mean(akt_err)
    }
  }
  ERR <- colMeans(ERR)
  # return a smallest local minimum of ERR
  return(min(localMaxima(-ERR)))
}

CV2 <- function(X,Folds=10,maxR=7,maxiter=10){
# performs CV for the best R (between 1 and maxR) like for bandwidth selection in KDE
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Nnew <- (N %/% Folds)*Folds
  Ind <- matrix(sample(1:Nnew), nrow=Folds)
  ERR <- array(0,c(Folds,maxR))
  Norms <- rep(0,maxR)
  for(fold in 1:Folds){
    Xtrain <- X[-Ind[fold,],,]
    Xtest <- X[Ind[fold,],,]
    Ntest <- dim(Xtest)[1]
    
    Res <- separable_expansion(Xtrain,maxR,maxiter)
    for(r in 1:maxR){
      Norms[r] <- sum(Res$sigma[1:r]^2)
      akt_err <- rep(0,Ntest)
      akt_Res <- Res
      akt_Res$A <- Res$A[1:r,,,drop=F]
      akt_Res$B <- Res$B[1:r,,,drop=F]
      akt_Res$sigma <- Res$sigma[1:r]
      for(n in 1:Ntest){
        akt_err[n] <- sum(apply_lhs(akt_Res,Xtest[n,,])*Xtest[n,,])
      }
      ERR[fold,r] <- mean(akt_err)
    }
  }
  ERR <- colMeans(ERR)
  CVscores <- Norms - 2*ERR
  # return a smallest local minimum of ERR
  return(min(localMaxima(-CVscores)))
}

###############
### Inverse ###
###############

eigvals <- function(Est,tol=1e-10,maxiter=10000){
# calculates the largest and the smallest eigenvalue of the estimated covariance
  K1 <- dim(Est$A)[2]
  K2 <- dim(Est$B)[2]
  V <- matrix(runif(K1*K2),ncol=K2) # initial guess
  V <- V/frobenius(V)
  eps <- 1
  iter <- 0
  while(eps > tol && iter < maxiter){
    iter <- iter+1
    V_old <- V
    V <- apply_lhs(Est,V)
    V <- V/frobenius(V)
    eps <- frobenius(V - V_old)
  }
  lambda_max <- sum(V*apply_lhs(Est,V))
  
  V <- matrix(runif(K1*K2),ncol=K2) # initial guess
  V <- V/frobenius(V)
  eps <- 1
  while(eps > tol && iter < maxiter){
    iter <- iter+1
    V_old <- V
    V <- apply_lhs(Est,V)
    V <- lambda_max*V_old - V
    V <- V/frobenius(V)
    eps <- frobenius(V - V_old)
  }
  lambda_min <- sum(V*apply_lhs(Est,V))
  return(list(min=lambda_min, max=lambda_max))
}

apply_lhs <- function(Est,Y){
# OUTPUT: Z = 'Est(Y)' = A1 Y B1 + ... + AR Y BR
  R <- length(Est$sigma)
  K1 <- dim(Est$A)[2]
  K2 <- dim(Est$B)[3]
  Z <- array(0,c(K1,K2))
  for(r in 1:R){
    Z <- Z + Est$sigma[r] * ( Est$A[r,,] %*% Y %*% Est$B[r,,] )
  }
  return(Z)
}

pcg <- function(Est,Y,Ig=NULL,tol=1e-10){
# PCG method for solving A1 X B1 + ... + AR X BR = Y with preconditioner (A1 x B1)^-1
# Est - a list giving the left-hand side:
#  $A - array[R,K1,K1], e.g. Est$A[r,,] gives matrix Ar
#  $B - array[R,K2,K2], e.g. Est$B[r,,] gives matrix Br
#  $sigma - vector[R] of singular values (A's nad B's standardized to norm 1)
# Y - the right-hand side matrix of size K1 x K2
# Ig - the initial guess (matrix of the same size a Y)
# tol - the tolerance
# OUTPUT: the solution as a matrix X
  K1 <- dim(Y)[1]
  K2 <- dim(Y)[2]
  mmax = 2000; # the maximal number of iterations
  if(length(Ig)==0) Ig <- matrix(rnorm(K1*K2),ncol=K2)
  U <- Ig
  R <- Y-apply_lhs(Est,U) # residual here, not Kronecker rank
  E <- rep(0,mmax)
  E[1] <- sqrt(sum(R^2))
  # cat('Initial residual',E[1],'\n')
  iter = 1
  t1 = 1.0
  D = array(0,c(K1,K2))
  A1_inv <- solve(Est$A[1,,])
  B1_inv <- solve(Est$B[1,,])
  while(iter < mmax && E[iter]/E[1] > tol){
    Z = ( A1_inv %*% R %*% B1_inv ) / Est$sigma[1]
    t1old = t1
    t1 = sum(Z*R) 
    beta = t1/t1old
    D = Z+beta*D
    S = apply_lhs(Est,D)
    suma = sum(D*S)
    tau = t1/suma
    U = U+tau*D
    R = R-tau*S
    iter = iter+1
    E[iter] = sqrt(sum(R^2))
    # if(iter %% 10 == 0) cat('Step', iter, 'relative residual', E[iter]/E[1],'\n')
  }
  return(list(U=U,iter=iter))
}

pcg2 <- function(Est,Y,eps,Ig=NULL,tol=1e-6){
  # PCG method for solving A1 X B1 + ... + AR X BR + eps X = Y with preconditioner (A1 x B1 + eps I)^-1
  # Est - a list giving the left-hand side:
  #  $A - array[R,K1,K1], e.g. Est$A[r,,] gives matrix Ar
  #  $B - array[R,K2,K2], e.g. Est$B[r,,] gives matrix Br
  #  $sigma - vector[R] of singular values (A's nad B's standardized to norm 1)
  # Y - the right-hand side matrix of size K1 x K2
  # eps - regularization for positive definiteness
  # Ig - the initial guess (matrix of the same size a Y)
  # tol - the tolerance
  # OUTPUT: the solution as a matrix X
  K1 <- dim(Y)[1]
  K2 <- dim(Y)[2]
  
  R <- length(Est$sigma)
  Est_new <- list(sigma=rep(0,R+1),A=array(0,c(R+1,K1,K1)), B=array(0,c(R+1,K2,K2)))
  Est_new$sigma[1:R] <- Est$sigma
  Est_new$A[1:R,,] <- Est$A
  Est_new$B[1:R,,] <- Est$B
  Est_new$sigma[R+1] <- eps
  Est_new$A[R+1,,] <- diag(K1)
  Est_new$B[R+1,,] <- diag(K2)
  
  mmax = 2000; # the maximal number of iterations
  if(length(Ig)==0) Ig <- matrix(rnorm(K1*K2),ncol=K2)
  U <- Ig
  R <- Y-apply_lhs(Est_new,U) # residual here, not Kronecker rank
  E <- rep(0,mmax)
  E[1] <- sqrt(sum(R^2))
  # cat('Initial residual',E[1],'\n')
  iter = 1
  t1 = 1.0
  D = array(0,c(K1,K2))
  EIG <- eigen(Est_new$B[1,,])
  UU <- EIG$vectors
  alpha1 <- EIG$values
  EIG <- eigen(Est_new$A[1,,])
  VV <- EIG$vectors
  alpha2 <- EIG$values
  HH <- ( (alpha2 %*% t(alpha1))*Est_new$sigma[1] + eps )^(-1)
  while(iter < mmax && E[iter]/E[1] > tol){    
    Z <- VV %*% ( HH * (t(VV) %*% R %*% UU) ) %*% t(UU)     # analytic solution to the Stein's equation
    t1old = t1
    t1 = sum(Z*R) 
    beta = t1/t1old
    D = Z+beta*D
    S = apply_lhs(Est_new,D)
    suma = sum(D*S)
    tau = t1/suma
    U = U+tau*D
    R = R-tau*S
    iter = iter+1
    E[iter] = sqrt(sum(R^2))
    # if(iter %% 10 == 0) cat('Step', iter, 'relative residual', E[iter]/E[1],'\n')
  }
  return(list(U=U,iter=iter))
}

#################
### utilities ###
#################

frobenius <- function(X){
  return(sqrt(sum(X^2)))
}

frobenius_PCA <- function(Res,A,B,R=NULL){
  R_true <- dim(A)[1]
  if(is.null(R)) R <- dim(Res$A)[1]
  K1 <- dim(Res$A)[2]
  K2 <- dim(Res$B)[2]
  C_hat <- Res$sigma[1]*outer(Res$A[1,,],Res$B[1,,])
  if(R>1){
    for(r in 2:R){
      C_hat <- C_hat + Res$sigma[r]*outer(Res$A[r,,],Res$B[r,,])
    }
  }
  C <- outer(A[1,,],B[1,,])
  if(R_true>1){
    for(r in 2:R_true){
      C <- C + outer(A[r,,],B[r,,])
    }
  }
  return(frobenius(C-C_hat)/frobenius(C))
  # return(sqrt(sum((C-C_hat)^2))/sqrt(sum(C^2)))
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

VanLoansPerm <- function(C_mat,K1,K2){
  return( tensor2matrix(aperm(matrix2tensor(C_mat,K1,K2),c(1,3,2,4))) )
}

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  return(y)
}

###################
### Simulations ###
###################

rnd_cov <- function(K=NULL,nf=NULL,sconst=NULL,decay=NULL,seed=NULL){
  # produces a random covariance with leading smoooth eigenfunctions from a certain list of prescribed functions
  # INPUT: K - grid size (size of the covariance matrix to be produced)
  #       nf - how many smooth eigenfunctions should be picked (has to be lower than 20)
  #   sconst - if large, eigenfunctions are really smooth, if low, they are perturbed (should be between 5 and 20 or so)
  #    decay - eigendecay costant, 2.72 corresponds to exp decay (should be between 1 and 2.72 or so)
  #     seed - seed for set.seed()
  # NOTE: to adjust smoothness in a random fashion, change the limits in the first 3 lines
  # OUTPUT: a K x K covariance matrix of frobenius norm 1
  if(is.null(nf)) nf <- sample(2:10,1)
  if(is.null(sconst)) sconst <- runif(1,10,25)
  if(is.null(decay)) decay <- runif(1,1.5,2.7)
  if(is.null(K)) K <- 50
  x <- seq(-(K-1)/2,(K-1)/2,1)/(K/2)
  # define candidate eigenfunctions on the grid from which you will pick later
  Cands <- array(0,c(K,20))
  Cands[,1] <- 1/sqrt(K)
  Cands[,2] <- x/frobenius(x)
  Cands[,3] <- x^2/frobenius(x^2)
  Cands[,4] <- x^3/frobenius(x^3)
  Cands[,5] <- x*(1-x)/frobenius(x*(1-x))
  Cands[,6] <- (1-x^2)/frobenius(1-x^2)
  Cands[,7] <- cos(pi*x)/frobenius(cos(pi*x))
  Cands[,8] <- sin(2*pi*x)/frobenius(sin(2*pi*x))
  Cands[,9] <- cos(2*pi*x)/frobenius(cos(2*pi*x))
  Cands[,10] <- sin(3*pi*x)/frobenius(sin(3*pi*x))
  Cands[,11] <- cos(3*pi*x)/frobenius(cos(3*pi*x))
  Cands[,12] <- sin(pi*x)/frobenius(sin(pi*x))
  Cands[,13] <- x*sin(pi*x)/frobenius(x*sin(pi*x))
  Cands[,14] <- x*cos(pi*x)/frobenius(x*cos(pi*x))
  Cands[,15] <- x*sin(2*pi*x)/frobenius(x*sin(2*pi*x))
  Cands[,16] <- x*cos(2*pi*x)/frobenius(x*cos(2*pi*x))
  y <- bs(x,degree=5)
  Cands[,17] <- y[,1]
  Cands[,18] <- y[,2]
  Cands[,19] <- y[,3]
  Cands[,20] <- y[,4]
  #
  C <- array(0,c(K,K))
  if(!is.null(seed)) set.seed(seed)
  picks <- sample(c(1:20),nf)
  C <- matrix(rnorm(K^2),ncol=K)
  for(j in 1:nf){
    C[,j] <- Cands[,picks[j]]
  }
  C[,1:nf] <- C[,1:nf]*sconst
  SVD <- svd(C)
  C <- SVD$u %*% diag(decay^(K:1)) %*% t(SVD$u)
  if(runif(1)>1/2) C <- apply(t(apply(C,2,rev)),2,rev)
  return(C/frobenius(C))
}

gneiting_cov <- function(t1,s1,t2,s2,beta,K1=NULL,K2=NULL){
# single point of the covariance introduced in Gneiting (2002) with parameters fixed as in Bagchi & Dette (2019)
# if K1 and K2 provided, t1,s1,t2,s2 can be integers, subsequently converted in equidistant grid on [0,1]
  if(length(K1)>0){
    t1 <- t1/(K1+1)
    t2 <- t2/(K1+1)
    s1 <- s1/(K2+1)
    s2 <- s2/(K2+1)
  }
  gamma <- 1 
  alpha <- 1
  sigma2 <- 1
  a <- 400 
  const <- 400 
  tau <- 1
  
  Cov <- const*abs(s1-s2)^(2*gamma)/(a*abs(t1-t2)^(2*alpha)+1)^(beta*gamma)
  Cov <- sigma2/(a*abs(t1-t2)^(2*alpha)+1)^tau * exp(-Cov)
  return(Cov)
}

get_gneiting_cov <- function(K1,K2,beta){
# returns Gneiting's covariance evaluated on the grid K1xK2 as a tensor
  C <- array(0,c(K1,K2,K1,K2))
  for(i in 1:K1){
    for(j in 1:K2){
      for(k in 1:K1){
        for(l in 1:K2){
          C[i,j,k,l] <- gneiting_cov(i,j,k,l,beta,K1,K2)
        }
      }
    }
  }
  return(C)
}

get_gneiting_symbol <- function(K1,K2,beta){
  # returns Gneiting's covariance evaluated on the grid K1xK2 as a tensor
  Gamma <- array(0,c(K1,K2))
  for(i in 1:K1){
    for(j in 1:K2){
      Gamma[i,j] <- gneiting_cov(i-1,j-1,0,0,beta,K1,K2)
    }
  }
  return(Gamma)
}

brownian_cov <- function(K){
  B <- array(0,c(K,K))
  for(i in 1:K){
    for(j in 1:K){
      B[i,j] <- min(i,j)
    }
  }
  return(B)
}
