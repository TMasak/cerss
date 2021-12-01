LDA <- function(mu0,mu1,Res,Xnew,psi,prior=T){
  a <- sum( (Xnew-mu0)*psi )^2
  b <- sum( (Xnew-mu1)*psi )^2
  if(prior){
    pi0 <- 45/121; pi1 <- 76/121; # this is what we have among our 121 patients
    const <- sum(psi*apply_lhs(Res,psi)) # should application here be regularized?
    return( I(a-b > 2*const*log(pi0/pi1)) )
  }else{
    return( I(a > b) )
  }
}

LDA_pca <- function(mu0,mu1,Res,Xnew,psi,prior=T){
  a <- sum( (Xnew-mu0)*psi )^2
  b <- sum( (Xnew-mu1)*psi )^2
  q <- length(Res$l)
  if(prior){
    pi0 <- 45/121; pi1 <- 76/121; # this is what we have among our 121 patients
    # apply_lhs
    res <- array(0,dim(psi))
    for(j in 1:q){
      res <- res + sum(Res$v[j,,]*psi)/Res$l[j]*Res$v[j,,]
    }
    const <- sum(psi*res) 
    return( I(a-b > 2*const*log(pi0/pi1)) )
  }else{
    return( I(a > b) )
  }
}

LDA_emp <- function(mu0,mu1,X,Xnew,psi,prior=T){
  a <- sum( (Xnew-mu0)*psi )^2
  b <- sum( (Xnew-mu1)*psi )^2
  if(prior){
    pi0 <- 45/121; pi1 <- 76/121; # this is what we have among our 121 patients
    const <- sum(psi*apply_lhs_emp(X,psi))
    return( I(a-b > 2*const*log(pi0/pi1)) )
  }else{
    return(I(a > b))
  }
}

PLS_psi <- function(mu0,mu1,Res,q){
  K1 <- dim(mu0)[1]
  K2 <- dim(mu0)[2]
  mu <- mu1 - mu0
  psi <- array(0,c(64,256))
  v <- e <- mu
  for(j in 1:q){
    h <- sum(v*e)/sum(v*apply_lhs(Res,v))
    psi <- psi + h*v
    e <- mu - apply_lhs(Res,psi)
    g <- sum(e*apply_lhs(Res,psi))/sum(v*apply_lhs(Res,v))
    v <- e + g*v
  }
  return(psi)
}

PLS_psi_emp <- function(mu0,mu1,X,q){
  K1 <- dim(mu0)[1]
  K2 <- dim(mu0)[2]
  mu <- mu1 - mu0
  psi <- array(0,c(64,256))
  v <- e <- mu
  for(j in 1:q){
    h <- sum(v*e)/sum(v*apply_lhs_emp(X,v))
    psi <- psi + h*v
    e <- mu - apply_lhs_emp(X,psi)
    g <- sum(e*apply_lhs_emp(X,psi))/sum(v*apply_lhs_emp(X,v))
    v <- e + g*v
  }
  return(psi)
}

apply_lhs_emp <- function(X,v){
  N <- dim(X)[1]
  res <- array(0,dim(v))
  for(n in 1:N){
    res <- res + sum(X[n,,]*v)*X[n,,]
  }
  return(res/N)
}

ridge_psi <- function(mu0,mu1,Res,q){
  K1 <- dim(mu0)[1]
  K2 <- dim(mu0)[2]
  mu <- mu1 - mu0
  # EIG <- eigvals(Res)
  # if(EIG$min < 0) eps <- -EIG$min + q else eps <- q
  eps <- q
  psi <- pcg(Res,mu,eps)
  psi <- psi$U
  return(psi)
}

power_iter <- function(X,q,tol=1e-6,maxiter=100){
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  V <- array(0,c(q,K1,K2))
  Iters <- lambda <- rep(0,q)
  for(j in 1:q){
    err <- iter <- 1
    v_old <- v <- array(runif(K1*K2),c(K1,K2))
    v <- v/base::norm(v,type="F")
    while(err > tol && iter < maxiter){
      v <- apply_lhs_emp(X,v)
      if(j > 1){ # deflate
        defval <- array(0,c(K1,K2))
        for(k in 1:(j-1)){
          defval <- defval + lambda[k]*sum(V[k,,]*v_old)*V[k,,]
        }
        v <- v - defval
      }
      v <- v/base::norm(v,type='F')
      err <- base::norm(v-v_old,type="F")
      v_old <- v
      iter <- iter + 1
    }
    V[j,,] <- v
    lambda[j] <- sum(v*apply_lhs_emp(X,v))
    Iters[j] <- iter
  }
  return(list(l=lambda,v=V,iter=Iters))
}

PCA_psi <- function(mu0,mu1,Res,q){
  K1 <- dim(mu0)[1]
  K2 <- dim(mu0)[2]
  mu <- mu1 - mu0
  psi <- array(0,c(K1,K2))
  for(j in 1:q){
    psi <- psi + sum(Res$v[j,,]*mu)/Res$l[j]*Res$v[j,,]
  }
  return(psi)
}

################################################################################

pcg <- function(Est,Y,eps,Ig=NULL,tol=1e-10){
  K1 <- dim(Y)[1]
  K2 <- dim(Y)[2]
  R <- length(Est$sigma)
  Est_new <- list(sigma=rep(0,R+1),A=array(0,c(R+1,K1,K1)),
                  B=array(0,c(R+1,K2,K2)))
  Est_new$sigma[1:R] <- Est$sigma
  Est_new$A[1:R,,] <- Est$A
  Est_new$B[1:R,,] <- Est$B
  Est_new$sigma[R+1] <- eps
  Est_new$A[R+1,,] <- diag(K1)
  Est_new$B[R+1,,] <- diag(K2)
  mmax = 2000; # the maximum number of iterations
  if(length(Ig)==0) Ig <- matrix(stats::rnorm(K1*K2),ncol=K2)
  U <- Ig
  R <- Y-apply_lhs(Est_new,U)
  E <- rep(0,mmax)
  E[1] <- sqrt(sum(R^2))
  if(!is.finite(E[1])){
    E[1] <- 1
    iter <- mmax
    U <- array(0,c(K1,K2))
  }
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
    # analytic solution to the Stein's equation
    Z <- VV %*% ( HH * (t(VV) %*% R %*% UU) ) %*% t(UU)
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
    if(!is.finite(E[iter])){
      E[iter] <- 0
      U <- array(0,c(K1,K2))
    }
  }
  return(list(U=U,iter=iter))
}

apply_lhs <- function(Est,Y){
  R <- length(Est$sigma)
  K1 <- dim(Est$A)[2]
  K2 <- dim(Est$B)[2]
  Z <- array(0,c(K1,K2))
  for(r in 1:R){
    Z <- Z + Est$sigma[r] * ( Est$A[r,,] %*% Y %*% Est$B[r,,] )
  }
  return(Z)
}

scd <- function(X, R=NULL, B=NULL, mu=NULL, predict=FALSE, maxiter=10,
                maxR=NULL, Folds=NULL, perc=c(1/4,1/4)){
  ### checking inputs
  # data set X
  if(length(dim(X)) != 3){
    stop(paste0("The data set has to form a 3-dimensional array."))
  }
  if(sum(I(dim(X) > 1)) < 3){
    stop(paste0("The data set has to form a 3-dimensional array with all
                dimensions being non-trivial."))
  }
  if(sum(!is.finite(X)) > 0){
    stop(paste0("Infinite or missing values in the data set."))
  }
  if(sum(!is.numeric(X)) > 0){
    stop(paste0("Non-numeric values in the data set."))
  }
  if(sum(I(apply(X, c(2,3), stats::var) < 10*.Machine$double.eps)) > 0){
    warning(paste0("Not enough variability in some locations."))
  }
  N <- dim(X)[1]
  Kmin <-  min(dim(X)[-1])
  # degree-of-separability R
  if(length(R) > 0){
    R <- as.integer(R)
    if(!is.integer(R) || R < 1){
      stop(paste0("The degree-of-separability R has to be a positive integer."))
    }
    if(R > N || R > Kmin)
    {
      stop(paste0("The chose degree-of-separability R is too high.
                     It needs to be R <= min(dim(X))."))
    }
  }
  # maximal degree for cross-validation maxR
  if(length(maxR)>0){
    maxR <- as.integer(maxR)
    if(!is.integer(maxR) || maxR < 1){
      stop(paste0("Them maximal degree-of-separability maxR has to be
                  a positive integer."))
    }
    if(maxR > N || maxR > Kmin)
    {
      stop(paste0("The chosen maximal degree-of-separability maxR is too high.
                  It needs to be maxR <= min(dim(X))"))
    }
  }
  # starting point - below, only when R is provided and no CV utilized
  # mean mu
  if(length(mu) > 0){
    if(dim(mu)[1] != dim(X)[2] || dim(mu)[2] != dim(X)[3]){
      warning(paste0("Dimensions of the mean mu do not correspond to the
                     data size. Discarding the provided mu."))
      mu <- NULL
    }
  }
  # prediction flag predict
  predict <- as.logical(predict)
  # maximum no. of iterations maxiter
  maxiter <- as.integer(maxiter)
  # number of splits for the cross-validation Folds
  if(length(Folds) > 0){
    Folds <- as.integer(Folds)
    Folds <- min(max(2,Folds),floor(N/2))
  }
  # percentages for prediction
  if(length(perc) > 1) perc <- perc[1:2] else perc <- c(perc,perc)
  if(perc[1] < 0 || perc[2] < 0){
    warning(paste0("Percentages for prediction need to be positive,
                   collapsing to the default."))
    perc <- c(1/4,1/4)
  }
  if(perc[1] >= Kmin || perc[2] >= Kmin){
    stop(paste0("Cannot understand the provided perc argument."))
  }
  if(perc[1] >= 1) perc[1] <- perc[1]/dim(X)[2]
  if(perc[2] >= 1) perc[2] <- perc[2]/dim(X)[3]
  ### set unused variables to their default
  if(length(mu)==0){ # use empirical mean to center the data
    mu <- apply(X, c(2,3), mean)
  }
  #
  if(length(maxR)==0) maxR <- 3
  if(length(Folds)==0) Folds <- min(10, floor(N/2))
  if(N < 2*Folds || Folds==1){
    stop(paste0("Not enough surfaces to do cross-validation, choose the
                bandwidth manually."))
  }
  ### cross-validation (CV)
  CV <- FALSE
  X <- sweep(X, c(2,3), mu)
  if(length(R)==0){
    CV <- TRUE
    if(predict){
      cvscores <- cvR_pred(X,Folds,maxR,maxiter,perc)
    }else{
      cvscores <- cvR(X,Folds,maxR,maxiter)
    }
    R <- min(localMaxima(-cvscores))
  }
  # check out the starting point B, if provided by the user
  if(length(B) > 0){
    if(length(dim(B)) != 3){
      if(length(dim(B)) == 2){
        Bb <- B
        B <- array(0,c(R,dim(B)))
        for(r in 1:R) B[r,,] <- Bb
      }else{
        warning(paste0("Cannot use the provided starting point B."))
        B <- NULL
      }
    }
    if(dim(B)[2] != dim(X)[3] || dim(B)[3] != dim(X)[3]){
      warning(paste0("Dimensions of the starting point B do not correspond to
                     the spatial covariance kernel of X. Discarding
                     the provided B."))
      B <- NULL
    }
  }
  ### Estimation
  if(CV){
    Res <- scd_est(X,R,maxiter)
    return(list(A=Res$A, B=Res$B, sigma=Res$sigma,
                cv=cvscores, R=R))
  }else return(scd_est(X,R,maxiter,B))
}

scd_est <- function(X,R=1,maxiter=10,B=NULL){
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
  for(r in 1:R){
    iter <- 0
    while(iter < maxiter){
      A[r,,] <- T1(X,A,B,sigma,r,N,K1,K2)
      A[r,,] <- A[r,,]/norm(A[r,,],type="F")
      B[r,,] <- T2(X,A,B,sigma,r,N,K1,K2)
      sigma[r] <- norm(B[r,,],type="F")
      B[r,,] <- B[r,,]/sigma[r]
      iter <- iter +1
    }
  }
  return(list(A=A,B=B,sigma=sigma))
}

T1 <- function (X,A,B,sigma,r,N,K1,K2){
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