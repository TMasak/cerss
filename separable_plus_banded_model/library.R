require(MASS)
require(orthopolynom)
# require(tensorA)

##################
### Estimation ###
##################

estimate_all <- function(Y, delta=0){
# Produces the separable plus banded estimator of the covariance
  Est <- estimate_separable(Y,delta=delta)
  B <- as.matrix(0)
  if(delta > 0){
    B <- as.matrix(toeplitz_average(Y,delta) - toeplitz_average_separable(Est$A1,Est$A2,delta))
  }
  return(list(A1 = Est$A1, A2 = Est$A2, B = B))
}

estimate_separable <- function(Y, delta=0, standardize=TRUE){
# INPUT:     Y - array of order N x K1 x K2 representing N observations of separable process observed on
#                a K1 x K2 grid
#        delta - discrete size of the corrupted band (delta=1 corresponds to corrupted diagonal)
# OUTPUT: 2D covariances A1 and A2
  N <- dim(Y)[1]
  K1 <- dim(Y)[2]
  K2 <- dim(Y)[3]
  W1 <- matrix(0,K1,K1) # weighting matrix, will have ones on the delta-off-diagonals (both off-diagonals,
  indx <- 1:(K1-delta)                                                     # making the results symmetric)
  W1[cbind(indx+delta,indx)] <- W1[cbind(indx,indx+delta)] <- 1
  W2 <- matrix(0,K2,K2) 
  indx <- 1:(K2-delta)
  W2[cbind(indx+delta,indx)] <- W2[cbind(indx,indx+delta)] <- 1
  Ycenter <- sweep(Y, c(2,3), apply(Y, c(2,3), mean)) # centering data
  
  A1 <- array(0, dim = c(K1,K1))
  A2 <- array(0, dim = c(K2,K2))
  for(n in 1:N){
    A1 <- A1 + ( Ycenter[n,,] %*% W2 %*% t(Ycenter[n,,]) )
    A2 <- A2 + ( t(Ycenter[n,,]) %*% W1 %*% Ycenter[n,,] )
  }
  A1 <- A1/N
  A2 <- A2/N
  A1 <- A1*sign(sum(diag(A1)))
  A2 <- A2*sign(sum(diag(A2))) # flip signs if shifted traces negative
  if(standardize){ 
    A1 <- A1/sqrt(abs(sum(A1*W1)))
    A2 <- A2/sqrt(abs(sum(A2*W2))) 
  }                        
  return(list(A1 = A1, A2 = A2))
}

### only works for square data at the moment
toeplitz_average <- function(Y, delta, unbiased=TRUE){
  # computes TopAvg(C_N), i.e. projection of empirical covariance operator (of Y) using 2D FFT
  # unbiased control the standardization, either ./(K-t)/(K-s) for TRUE or /K^2 for FALSE (this is biased)
  # OUTPUT - band - matrix of size delta x delta containing all information needed to
  #                 reconstruct b(,,,) into form of a matrix by function expand_b()
  N <- dim(Y)[1]  
  K <- dim(Y)[2]
  band <- array(0, dim = c(K, K))
  # Y <- sweep(Y, c(2,3), apply(Y, c(2,3), mean))
  for(n in 1:N){
    IFTFT <- fft(abs(fft(Y[n,,]))^2, inverse = TRUE)/K^2
    band <- band + Re(IFTFT)
  }
  for(t in 0:(delta-1)){
    for(s in 0:(delta-1)){
      if(unbiased) band[t+1,s+1] <- band[t+1,s+1]/(K-t)/(K-s) 
      else         band[t+1,s+1] <- band[t+1,s+1]/K^2
    }
  }
  return(band[1:delta,1:delta]/N)
}

toeplitz_average_separable <- function(L1, L2, delta, unbiased=TRUE){
# computes the [1:delta,1:delta] submatrix of TopAvg(A1 x A2)
# unbiased controls the standardization, either ./(K-t)/(K-s) for TRUE or /K^2 for FALSE (this is biased)
# OUTPUT - band - matrix of size delta x delta containing all information needed to
#                 reconstruct b(,,,) into a tensor by function expand_b()
  K <- dim(L1)[1]
  band <- array(0, dim = c(delta, delta))
  for(t in 0:(delta-1)){
    ind <- 1:(K-t)
    alpha <- sum(L1[cbind(ind+t,ind)])
    for(s in 0:(delta-1)){
      ind <- 1:(K-s)
      beta <- sum(L2[cbind(ind+s,ind)])
      if(unbiased) band[t+1,s+1] <- alpha*beta/(K-t)/(K-s) 
      else         band[t+1,s+1] <- alpha*beta/K^2
    }
  }
  return(band)
}

##################################
### ALS approach to estimation ###
##################################

PIP1 <- function(Ycenter,W1){
  N <- dim(Ycenter)[1]
  K <- dim(Ycenter)[2]
  covariance_2 <- array(0, dim = c(K,K))
  for(n in 1:N){
    covariance_2 <- covariance_2 + ( t(Ycenter[n,,]) %*% W1 %*% Ycenter[n,,] )
  }
  covariance_2 <- covariance_2/N
  return(covariance_2)
}

PIP2 <- function(Ycenter,W2){
  N <- dim(Ycenter)[1]
  K <- dim(Ycenter)[2]
  covariance_1 <- array(0, dim = c(K,K))
  for(n in 1:N){
    covariance_1 <- covariance_1 + ( Ycenter[n,,] %*% W2 %*% t(Ycenter[n,,]) )
  }
  covariance_1 <- covariance_1/N
  return(covariance_1)
}

PIP1_toeplitz <- function(band,W1){
  K <- dim(W1)[1]
  delta <- dim(band)[1]
  covariance_2 <- array(0,c(K,K))
  w <- rep(0,delta)
  for(i in 0:(delta-1)){
    ind <- 1:(K-i)
    w[i+1] <- sum(W1[cbind(ind+i,ind)])/(K-i)
  }
  for(j in 0:(delta-1)){
    ind <- 1:(K-j)
    covariance_2[cbind(ind+j,ind)] <- covariance_2[cbind(ind,ind+j)] <- sum(band[j+1,]*w)
  }
  return(covariance_2)
}

PIP2_toeplitz <- function(band,W2){
  K <- dim(W2)[1]
  delta <- dim(band)[1]
  covariance_1 <- array(0,c(K,K))
  w <- rep(0,delta)
  for(i in 0:(delta-1)){
    ind <- 1:(K-i)
    w[i+1] <- sum(W2[cbind(ind+i,ind)])/(K-i)
  }
  for(j in 0:(delta-1)){
    ind <- 1:(K-j)
    covariance_1[cbind(ind+j,ind)] <- covariance_1[cbind(ind,ind+j)] <- sum(band[,j+1]*w)
  }
  return(covariance_1)
}

ALS_estimate <- function(Y, L2=NULL, B=NULL, maxiter=100, tol=10^-6){
  N <- dim(Y)[1]
  K <- dim(Y)[2]
  if(length(L2)==0){ L2 <- array(1,c(K,K))/K }
  if(length(B)==0){ B <- array(0,c(K,K)) }
  L1 <- array(1,c(K,K))/K
  L1all <- array(0,c(maxiter,K,K))
  L2all <- array(0,c(maxiter,K,K))
  Ball <- array(0,c(maxiter,K,K))
  iter <- 0
  eps <- base::norm(L2,type="F")
  Ycenter <- sweep(Y, c(2,3), apply(Y, c(2,3), mean))
  while(eps > tol && iter < maxiter){
    iter <- iter + 1
    L1_old <- L1; L2_old <- L2; B_old <- B
    L1 <- PIP2(Ycenter, L2) - PIP2_toeplitz(B,L2)
    L1all[iter,,] <- L1 <- L1/base::norm(L1,type="F")
    L2all[iter,,] <- L2 <- PIP1(Ycenter, L1) - PIP1_toeplitz(B,L1)
    # standardize
    # L2all[iter,,] <- L2 <- L2/base::norm(L2,type="F")
    # pom <- PIP2(Ycenter, L2) - PIP2_toeplitz(B,L2)
    # pom <- sum(L1*pom)
    # L1all[iter,,] <- L1 <- L1*pom
    # L2all[iter,,] <- L2 <- L2*pom
    
    Binit <- toeplitz_average(Ycenter,K,FALSE) - toeplitz_average_separable(L1,L2,K,FALSE)
    ### Should TopAvg be called on Y or Ycenter ?
    # Binit <- estimate_stationary_band(Y,floor(K/2),N,K,L1,L2)
    # ? for some reason we have to zero-out B to cancel out DFT symmetry ?
    B <- array(0,c(K,K))
    B[1:floor(K/2),1:floor(K/2)] <- Binit[1:floor(K/2),1:floor(K/2)]
    Ball[iter,,] <- B
    # stopping criteria
    eps <- max( base::norm(L1 - L1_old, type="F")/base::norm(L1_old,type="F") ,
                base::norm(L2 - L2_old, type="F")/base::norm(L2_old,type="F") ,
                base::norm(B - B_old, type="F")/base::norm(B_old,type="F") )
  }
  return(list(A1 = L1all[1:iter,,], A2 = L2all[1:iter,,], B = Ball[1:iter,,], iter = iter))
}

ALS_separable <- function(Y, L2=NULL, maxiter=100, tol=10^-9){
  N <- dim(Y)[1]
  K <- dim(Y)[2]
  if(length(L2)==0){ L2 <- array(1,c(K,K)) }
  L1 <- array(1,c(K,K))
  L1all <- array(0,c(maxiter,K,K))
  L2all <- array(0,c(maxiter,K,K))
  iter <- 0
  eps <- base::norm(L2,type="F")
  Ycenter <- sweep(Y, c(2,3), apply(Y, c(2,3), mean))
  while(eps > tol && iter < maxiter){
    iter <- iter + 1
    L1_old <- L1; L2_old <- L2
    L1 <- PIP2(Ycenter, L2)
    L1all[iter,,] <- L1 <- L1/base::norm(L1,type="F") # fix norm to 1
    L2all[iter,,] <- L2 <- PIP1(Ycenter, L1)
    # stopping criteria
    eps <- max( base::norm(L1 - L1_old, type="F")/base::norm(L1_old,type="F") ,
                base::norm(L2 - L2_old, type="F")/base::norm(L2_old,type="F") )
  }
  return(list(L1 = L1all[1:iter,,], L2 = L2all[1:iter,,], iter = iter ))
}

#######################################
### ADI solution to inverse problem ###
#######################################

ADI <- function(A1, A2, band, y, theta=1e-5, rho=0, adapt=T, maxiter=200, tol=10^-9){
# Solution to the inverse problem involving separable+stationary LHS
# INPUT:  A1, A2 - separable constituents
#           band - symbol of the stationary part (only its non-zero part)
#              y - RHS as a vector
#            rho - initial shift parameter
#          theta - additional regularization
#          adapt - if TRUE, rho changes between iterations
#        maxiter - maximum number of ADI iterations
#            tol - tolerance for stopping criterion
# OUTPUT:       x - the solution as a vector
#            iter - no of ADI iterations
#         PCGiter - vector of length iter, giving no of PCG iterations inside every sinlge ADI iteration
  K <- dim(A1)[1]
  iter <- 0
  eps <- 1
  x <- runif(K^2)
  PCGiter <- rep(0,maxiter)
  # pre-calculate eigen-decompositions
  A <- A1; A1 <- A2; A2 <- A # kronecker product definition's artifact
  EIG <- eigen(A1)
  U <- EIG$vectors
  alpha1 <- EIG$values
  EIG <- eigen(A2)
  V <- EIG$vectors
  alpha2 <- EIG$values
  while(iter < maxiter && eps > tol){
    x_old <- x
    iter <- iter + 1
    # calculate x^{(k+1/2)}
    tev <- Re(fft(to_book_format(band,K)))
    y_temp <- y - tx(tev,x_old) + rho*x_old            # RHS for the first half of the iteration
    Y <- matrix(y_temp, ncol=K)                        #
    H <- ( (alpha2 %*% t(alpha1)) + rho + theta )^(-1) #
    X <- V %*% ( H * (t(V) %*% Y %*% U) ) %*% t(U)     # analytic solution to the Stein's equation
    x_half <- c(X)
    # calculate x^{(k+1)}
    y_temp <- y - c( A2 %*% X %*% A1 ) + rho*x_half
    
    band_reg <- band
    band_reg[1,1] <- band_reg[1,1] + rho + theta
    ev <- circ_eig(to_book_format(band_reg,K))
    tev <- Re(fft(to_book_format(band_reg,K)))
    PCG <- pcg(tev,y_temp,runif(K^2),ev,tol)
    x <- PCG$u
    PCGiter[iter] <- PCG$iter
    # x <- Matrix::solve(B+(rho+theta)*diag(K^2), y_temp)
    
    eps <- sum((x-x_old)^2)/sum(x_old^2)
    if(adapt){
      rho <- min(rho, eps)
    }
  }
  return(list(x=x, iter=iter, PCGiter=PCGiter[1:iter]))
}

to_book_format <- function(band,K){
# transforms our stationary banded format of B into the format suitable for 2D FFT
  Bnew <- array(0, c(2*K,2*K))
  delta <- dim(band)[1]
  Bnew[1:delta,1:delta] <- band
  Bnew[(2*K):(K+2),1:K] <- Bnew[2:K,1:K]
  Bnew[1:K,(2*K):(K+2)] <- Bnew[1:K,2:K]
  Bnew[(2*K):(K+2),(2*K):(K+2)] <- Bnew[2:K,2:K]
  return(Bnew)
}

cg <- function(tev,b,ig,tol=1e-7){
# standard CG method for solving B x = y without any preconditioner
# INPUT: tev - eigenvalues of B (fully determine B)
#          b - the right-hand side vector
#         ig - the initial guess for the CG method
#        tol - tolerance for stopping criterion
# OUTPUT: solution 'u' as a vector, no of iterations 'iter'
  mmax = 2000;                     # the maximal number of iterations
  u = ig;
  r = b-tx(tev,u)
  e <- rep(0,mmax)
  e[1] = sqrt(sum(r^2))
  # cat('Initial residual',e[1],'\n')
  iter <- 1
  t1 <- 1.0
  d <- rep(0, length(ig)) 
  while(iter < mmax && e[iter]/e[1] > tol){
    z =r
    t1old = t1
    t1 = sum(z*r)
    beta = t1/t1old
    d = z+beta*d
    s = tx(tev,d)
    suma = sum(d*s)
    tau = t1/suma
    u = u+tau*d;
    r = r-tau*s;
    iter = iter+1;
    e[iter] = sqrt(sum(r^2))
    # cat('Step', iter, 'relative residual', e[iter]/e[1],'\n')
  }
  return(list(u=u, iter=iter))
}

tx <- function(tev,v){
# fast multiplication B v
# tev - eigenvalues of B (fully determine B)
# OUTPUT: the product given in `y` as a vector
  n1 <- dim(tev)[1]
  m1 <-dim(tev)[2]
  n <- n1/2
  m <- m1/2
  v <- matrix(v, ncol=m)
  ev <- array(0,c(n1,m1))
  ev[1:n,1:m] <- v
  y <- fft(ev)
  y <- tev*y
  y <- fft(y, inverse = T)/m1/n1
  y <- y[1:n,1:m]
  return(c(Re(y)))
}

circ_eig <- function(t){
# Computes eigenvalues of blocks of the preconditioner (denoted c^{(1)}_F in Chan & Jin, 2007)
# INPUT: the matrix generated by to_book_format()
# OUTPUT: n-by-(2m) matrix, each column consists of eigenvalues of a circulant block.
  n <- dim(t)[1]
  m <- dim(t)[2]
  n = n/2
  ev = array(0,c(n,m))
  v = array(0,c(n,m))
  v[1,] = t[1,];
  for (i in 2:n) v[i,] <- ((n-(i-1))*t[i,]+(i-1)*t[n+i,])/n;
  ev <- mvfft(v)  # FFT applied column by column (not multidimensional FFT)
  return(Re(ev))
}

psolve <- function(ev,d){
# Solves the preconditioning system Cy=d (changes back the variables).
# INPUT: ev - eigenvalues of the preconditioner, generated by circ_eig (fully determine C)
#         d - the right-hand side vector
# OUTPUT - the solution `y`
  n <- dim(ev)[1]
  m <- dim(ev)[2]
  m <- m/2
  rex = matrix(d,ncol=m)
  rex = mvfft(rex)      # FFT column by column
  for(i in 1:n){
    A = toeplitz(ev[i,1:m])
    rex[i,] = t(solve( A, rex[i,] ))
    # We may solve the Toeplitz systems by the PCG methods or
    # by fast direct methods
  }
  rex = mvfft(rex, inverse = T)/n # inverse FFT column by column
  return(c(Re(rex)))
}

pcg <- function(tev,b,ig,ev,tol=1e-7){
# PCG method for solving B x = b with the circulant preconditioner of Chan & Yin (2007)
# INPUT: tev - eigenvalues of B, fully specifying B
#          b - the right-hand side vector
#         ig - the initial guess
#        tol - the tolerance for stopping criterion
#         ev - eigenvalues of the circulant preconditioner
# OUTPUT: the solution `u`, no of iterations needed `iter`
  mmax = 2000; # the maximal number of iterations
  u <- ig
  r <- b-tx(tev,u)
  e <- rep(0,mmax)
  e[1] <- sqrt(sum(r^2))
  # cat('Initial residual',e[1],'\n')
  iter = 1
  t1 = 1.0
  d = rep(0, length(ig))
  while(iter < mmax && e[iter]/e[1] > tol){
    z = psolve(ev,r)
    t1old = t1
    t1 = sum(z*r)
    beta = t1/t1old
    d = z+beta*d
    s = tx(tev,d)
    suma = sum(d*s)
    tau = t1/suma
    u = u+tau*d
    r = r-tau*s
    iter = iter+1
    e[iter] = sqrt(sum(r^2))
    # cat('Step', iter, 'relative residual', e[iter]/e[1],'\n')
  }
  # if (iter == mmax) cat('Maximum iterations reached')
  return(list(u=u, iter=iter))
}

ADI_var <- function(A1, A2, B_var, y, theta=0, rho=0, adapt=F, maxiter=200, tol=10^-6){
# just as ADI is for separable+stationary LHS, this is for separable+diagonal LHS, i.e. for
# (separable + heteroscedastic variance) model
  K1 <- dim(A1)[1]
  K2 <- dim(A2)[1]
  iter <- 0
  eps <- 1
  x <- runif(K1*K2)
  # pre-calculate eigen-decompositions
  A <- A1; A1 <- A2; A2 <- A # kronecker product definition's artifact
  EIG <- eigen(A1)
  U <- EIG$vectors
  alpha1 <- EIG$values
  EIG <- eigen(A2)
  V <- EIG$vectors
  alpha2 <- EIG$values
  while(iter < maxiter && eps > tol){
    x_old <- x
    iter <- iter + 1
    # calculate x^{(k+1/2)}
    y_temp <- y - B_var*x_old + rho*x_old               # RHS for the first half of the iteration
    Y <- matrix(y_temp, ncol=K2)                        #
    H <- ( (alpha2 %*% t(alpha1)) + rho + theta )^(-1)  #
    X_half <- V %*% ( H * (t(V) %*% Y %*% U) ) %*% t(U) # analytic solution to the Stein's equation
    # calculate x^{(k+1)}
    y_temp <- y - c( A2 %*% X_half %*% A1 ) + rho*c(X_half)
    Y <- matrix(y_temp, ncol=K2)
    B_akt <- (B_var + rho)^(-1)
    X <- B_akt*Y
    x <- c(X)
    
    eps <- sum((x-x_old)^2)/sum(x_old^2)
    if(adapt){
      rho <- min(rho, eps)
    }
  }
  return(list(x=x, iter=iter))
}

######################################
### prediction and delta selection ###
######################################

CV <- function(X,Folds=10,maxd=20,mind=1){
# performs CV for the best delta (between 1 and maxd) like for bandwidth selection in KDE
  N <- dim(X)[1]
  K <- K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Nnew <- (N %/% Folds)*Folds
  Ind <- matrix(sample(1:Nnew), nrow=Folds)
  Fit <- array(0,c(Folds,maxd-mind+1))
  Norms <- rep(0,maxd-mind+1)
  for(delta in mind:maxd){
    Chat <- estimate_all(X,delta)
    Norms[delta-mind+1] <- frobenius_all(Chat,relative=F)
  }
  for(fold in 1:Folds){
    # print(fold)
    Xtrain <- X[-Ind[fold,],,]
    Xtest <- X[Ind[fold,],,]
    Ntest <- dim(Xtest)[1]
    for(delta in mind:maxd){
      Chat <- estimate_all(X,delta)
      ### I calculated norms here at first
      akt_err <- rep(0,Ntest)
      for(n in 1:Ntest){
        tev <- Re(fft(to_book_format(Chat$B,K)))
        pom <- Chat$A1 %*% Xtest[n,,] %*% Chat$A2
        pom <- pom + tx(tev,c(Xtest[n,,]))
        akt_err[n] <- sum(pom*Xtest[n,,])
      }
      Fit[fold,delta-mind+1] <- mean(akt_err)
    }
  }
  Fit <- colMeans(Fit)
  
  RES <- array(0,c(2,maxd-mind+1)) # 5...{Fit, Norms, sepFit, sepNorms, statNorms}
  RES[1,] <- Fit
  RES[2,] <- Norms
  return(RES)
}

stability_selection <- function(Y, maxd=20, mind=1){
# delta choice based on stability of the solution
  K <- dim(Y)[2]
  errA <- rep(0,maxd-mind)
  errB <- rep(0,maxd-mind)
  err <- rep(0,maxd-mind)
  delta <- mind
  Cnew <- Cold <- list(A1=0,A2=0,B=0)
  Cold <- estimate_all(Y,delta)
  RES <- array(0,c(3,maxd-mind))
  for(delta in (mind+1):maxd){
    # print(delta)
    Cnew <- estimate_all(Y,delta)
    errA[delta-mind] <- frobenius_separable(Cold$A1,Cold$A2,Cnew$A1,Cnew$A2)
    errB[delta-mind] <- frobenius_stationary(K,Cold$B,Cnew$B)
    err[delta-mind] <- frobenius_all(Cnew,Cold,F)
    Cold <- Cnew
  }
  RES[1,] <- errA
  RES[2,] <- errB
  RES[3,] <- err
  return(RES)
}

CV_nonstationary <- function(X,Folds=10,maxd=20,mind=1){
  # performs CV for the best delta (between 1 and maxd) like for bandwidth selection in KDE
  X <- sweep(X, c(2,3), apply(X, c(2,3), mean))
  N <- dim(X)[1]
  K <- K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Nnew <- (N %/% Folds)*Folds
  Ind <- matrix(sample(1:Nnew), nrow=Folds)
  Fit <- array(0,c(Folds,maxd-mind+1))
  Norms <- rep(0,maxd-mind+1)
  for(delta in mind:maxd){
    Chat <- estimate_separable(X,delta)
    C <- calculate_covariance_tensor(X) - aperm(outer(Chat$A1, Chat$A2),c(1,3,2,4))
    for(i in 1:dim(C)[1]){
      for(j in 1:dim(C)[2]){
        for(k in 1:dim(C)[3]){
          for(l in 1:dim(C)[4]){
            if(abs(i-j) >= delta || abs(k-l) >= delta) C[i,j,k,l] <- 0
          }
        }
      }
    }
    Norms[delta-mind+1] <- frobenius(C + aperm(outer(Chat$A1, Chat$A2),c(1,3,2,4)))
  }
  for(fold in 1:Folds){
    print(fold)
    Xtrain <- X[-Ind[fold,],,]
    Xtest <- X[Ind[fold,],,]
    Ntest <- dim(Xtest)[1]
    for(delta in mind:maxd){
      Chat <- estimate_separable(Xtrain,delta)
      C <- calculate_covariance_tensor(Xtrain) - aperm(outer(Chat$A1, Chat$A2),c(1,3,2,4))
      for(i in 1:dim(C)[1]){
        for(j in 1:dim(C)[2]){
          for(k in 1:dim(C)[3]){
            for(l in 1:dim(C)[4]){
              if(abs(i-j) >= delta || abs(k-l) >= delta) C[i,j,k,l] <- 0
            }
          }
        }
      }
      ### I calculated norms here at first
      akt_err <- rep(0,Ntest)
      for(n in 1:Ntest){
        pom <- Chat$A1 %*% Xtest[n,,] %*% Chat$A2
        pom <- c(pom) + to_matrix(C) %*% c(Xtest[n,,])
        akt_err[n] <- sum(pom*c(Xtest[n,,]))
      }
      Fit[fold,delta-mind+1] <- mean(akt_err)
    }
  }
  Fit <- colMeans(Fit)
  
  RES <- array(0,c(2,maxd-mind+1)) # 5...{Fit, Norms, sepFit, sepNorms, statNorms}
  RES[1,] <- Fit
  RES[2,] <- Norms
  return(RES)
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

diagonal_band <- function(X,A1,A2){
# calculates the diagonal of B without the assumption of stationarity
  K1 <- dim(A1)[1]
  K2 <- dim(A2)[1]
  B_var <- array(0,c(K1,K2))
  for(k1 in 1:K1){
    for(k2 in 1:K2){
      B_var[k1,k2] <- var(X[,k1,k2]) - A1[k1,k1]*A2[k2,k2]
    }
  }
  return(B_var)
}

#################
### utilities ###
#################

frobenius <- function(X){
# Calculates the Frobenius norm of an array (number, vector, matrix or a tensor)
  return(sqrt(sum(X^2)))
}

frobenius_all <- function(Chat, C=NULL, relative=T){
# calculate frobenius(outer(A1,B1) + Banded1 - outer(A2,B2) - Banded2) without the necessity store
# the large outer products or frobenius(outer(A1,B1)) + Banded1 if only two matrices are given  
# Banded stationary matrices only given by their bands
  K <- dim(Chat$A1)[1]
  FROB <- 0
  FROB_down <- 0
  delta1 <- dim(Chat$B)[1]
  if(is.null(C)){
    A1 <- A2 <- matrix(0, K, K)
    B <- matrix(0, delta1, delta1)
    C <- list(A1=A1, A2=A2, B=B)
  }
  delta2 <- dim(C$B)[1]
  Band1 <- matrix(0, K, K)
  Band2 <- matrix(0, K, K)
  Band1[1:delta1,1:delta1] <- Chat$B
  Band2[1:delta2,1:delta2] <- C$B
  for(i in 0:(K-1)){
    ind <- 1:(K-i)
    a1 <- Chat$A1[cbind(ind+i,ind)]
    a2 <- C$A1[cbind(ind+i,ind)]
    for(j in 0:(K-1)){
      ind <- 1:(K-j)
      b1 <- Chat$A2[cbind(ind+j,ind)]
      b2 <- C$A2[cbind(ind+j,ind)]
      beta <- Band1[i+1,j+1] - Band2[i+1,j+1]
      times <- 2
      if(i==0 && j==0) times <- 1
      if(i>0 && j>0) times <- 4
      FROB <- FROB + times*frobenius_piece(a1,b1,a2,b2,beta)
      FROB_down <- FROB_down + times*frobenius_piece(a2,b2,NULL,NULL,Band2[i+1,j+1])
    }
  }
  if(relative) FROB <- sqrt(FROB/FROB_down)
  else FROB <- sqrt(FROB)
  return(FROB)
}

frobenius_piece <- function(a1,b1,a2,b2,beta){
  if(length(a2)==0) return(frobenius(a1 %*% t(b1) + beta)^2)
  else return(frobenius(a1 %*% t(b1) - a2 %*% t(b2) + beta)^2)
}

frobenius_piece_fast <- function(a1,b1,a2,b2,beta){
# Calculates the frobenius norm^2 of a single K x K submatrix of mat(C) in R^{K^2 x K^2} in O(K) flops
# using Gram-Schmidt orthogonalization. It only makes sense to use this if K is very large.
  K1 <- length(a1)
  K2 <- length(b1)
  if(is.null(a2)){
    a2 <- rep(0,K1)
    b2 <- rep(0,K2)
  }
  M <- array(0,c(K1,3))
  M[,1] <- a1
  M[,2] <- a2
  M[,3] <- sqrt(beta)
  pom <- qr(M)
  S1 <- pom$qr[1:3,1:3]
  Ind <- lower.tri(S1)
  S1[Ind] <- 0
  
  M <- array(0,c(K2,3))
  M[,1] <- b1
  M[,2] <- -b2
  M[,3] <- sqrt(beta)
  pom <- qr(M)
  S2 <- pom$qr[1:3,1:3]
  Ind <- lower.tri(S2)
  S2[Ind] <- 0
  return(frobenius( S1 %*% t(S2) )^2)
}

frobenius_separable <- function(A1,B1,A2=NULL,B2=NULL){
  # calculate frobenius(outer(A1,B1) - outer(A2,B2)) without the necessity store the large outer products
  # or        frobenius(outer(A1,B1)) if only two matrices are given  
  n1 <- dim(A1)[1]
  n2 <- dim(A1)[2]
  FROB <- 0
  if(length(A2)==0){
    A2 <- matrix(0, n1, n2)
    B2 <- matrix(0, dim(B1)[1], dim(B1)[2])
  }
  for(i in 1:n1){
    for(j in 1:n2){
      FROB <- FROB + sum((A1[i,j]*B1 - A2[i,j]*B2)^2)
    }
  }
  return(sqrt(sum(FROB)))
}

frobenius_stationary <- function(K,A,B=NULL){
  # calculates frobenius error of the banded part
  delta <- dim(A)[1]
  FROB <- 0
  if(length(B)==0){
    B <- matrix(0, delta, delta)
  }
  delta2 <- dim(B)[1]
  if(delta2 < delta){
    C <- array(0,c(delta,delta))
    C[1:delta2,1:delta2] <- B
    B <- C
  }
  if(delta2 > delta){
    C <- array(0,c(delta2,delta2))
    C[1:delta,1:delta] <- A
    A <- C
    delta <- delta2
  }
  W <- array(0,c(delta,delta))
  if(delta>1){
    for(a in 0:(delta-1)){
      for(b in 0:(delta-1)){
        if(a==0 || b==0){ W[a+1,b+1] <- 2*(K-a)*(K-b) }
        else{ W[a+1,b+1] <- 4*(K-a)*(K-b) } 
      }
    }  
  }
  W[1,1] <- K^2
  sqrt(sum(W*(A-B)^2))
}

frobenius_empirical <- function(Y, A, B, band){
  Y <- sweep(Y, c(2,3), apply(Y, c(2,3), mean))
  N <- dim(Y)[1]
  K <- dim(A)[1]
  delta <- dim(band)[1]
  FROB <- 0
  FROB_down <- 0
  Band <- matrix(0, K, K)
  Band[1:delta,1:delta] <- band
  for(i in 1:K){
    for(j in 1:K){
      for(k in 1:K){
        for(l in 1:K){
          C_entry <- sum(Y[,i,k]*Y[,j,l])/N # C_N[i,j,k,l]
          FROB <- FROB + ( A[i,j]*B[k,l] + Band[1+abs(i-j),1+abs(k-l)] - C_entry )^2
          FROB_down <- FROB_down + ( A[i,j]*B[k,l] + Band[1+abs(i-j),1+abs(k-l)] )^2
        }
      }
    }
  }
  return(sqrt(sum(FROB)/sum(FROB_down)))
}

rootmatrix <- function(x){
# Calculates "square-root" of a given matrix, taken from the elasticnet package. 
  x.eigen <- eigen(x)
  d <- x.eigen$values
  d <- (d + abs(d))/2
  v <- x.eigen$vectors
  return(v %*% diag(sqrt(d)) %*% t(v))
}

to_matrix <- function(X){
# change a given 4D covariance array X to a covariance matrix matrix
  d <- dim(X)
  v <- c(X)
  M <- matrix(v, ncol=d[1]*d[2])
  return(M)
}

to_array <- function(X,d){
# inverse of to_matrix, here we need to know the dimensions of the array d=c(d[1],d[2])
  v <- c(X)
  A <- array(v, dim=c(d[1],d[2],d[1],d[2]))
  return(A)
}

expand_b <- function(band,delta,K){
# reconstructs the output of estimate_stationary_band() into the tensor corresponding to the banded part
  X <- array(0, c(K,K,K,K))
  indx <- rep(1:K,K)
  indy <- rep(1:K,each=K)
  for(a in 0:(delta-1)){
    for(b in 0:(delta-1)){
      indx <- rep(1:(K-a),K-b)
      indy <- rep(1:(K-b),each=K-a)
      X[cbind(indx+a,indy+b,indx,indy)] <- band[a+1,b+1] 
      X[cbind(indx,indy,indx+a,indy+b)] <- band[a+1,b+1] 
      X[cbind(indx+a,indy,indx,indy+b)] <- band[a+1,b+1] 
      X[cbind(indx,indy+b,indx+a,indy)] <- band[a+1,b+1] 
    }
  }
  return(X)
}

localMaxima <- function(x) {
# gives indices of local maxima of vector x
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  return(y)
}

max_prominence2 <- function(x){
# Returns index of the most prominent local minimum of a vector, where prominence is defined as
# height of the valley.
  xmin <- localMaxima(-x)
  m <- length(xmin)
  if(m==1) return(xmin)  # if there is only one local minimum, return that one
  prominence <- rep(0,m)
  for(j in 1:m){
    prominence[j] <- prominence[j] + max(x[xmin[j]:length(x)]) - x[xmin[j]]
    prominence[j] <- prominence[j] + max(x[1:xmin[j]]) - x[xmin[j]]
  }
  return(xmin[which.max(prominence)])
}

#################################
### Nearest Kronecker Product ###
#################################

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

###################
### Simulations ###
###################

shifted_legendre <- function(rank, K, N, normalized=TRUE){
# Generates centered observations with covariance having Legendre polynomials as eigenfunctions
# INPUT: rank - order of the highest polynomial (so together with the constant, true rank is n+1)!!!
#           K - grid size
#           N - no. of realizations
#  normalized - TRUE if polynomials should be normalized  
# OUTPUT: matrix with observations as columns
  rank <- rank-1
  polynomials <- slegendre.polynomials( rank, normalized )
  
  grid <- seq(from=0,to=1-1/K,by=1/K) + 1/(2*K)
  pol_val <- polynomial.values(polynomials, grid)
  pol_val <- matrix(unlist(pol_val),nrow=K)
  
  rank <- rank+1
  C <- matrix(rnorm(N*rank,0,1),nrow=rank)
  lb <- 0.5
  ub <- 1.2
  lambda <- sort(seq(lb,ub,length.out=rank),decreasing = TRUE)
  C <- C*lambda
  pol_val <- t(pol_val%*%C)
  pol_val <- t(t(pol_val)-apply(pol_val,2,mean))
  return(pol_val)
}

create_banded_datum <- function(mask,K){
# Creates a datum with stationary covariance B (given by `mask`).
# Fast way to do space-time averaging using FFT.
  delta <- dim(mask)[1]
  # mask <- mask/frobenius(mask)/K # to have trace fixed to one
  q <- array(0,c(K+delta-1,K+delta-1))
  q[1:delta,1:delta] <- mask
  q <- c(q)
  q <- c(q[1],rev(q[2:length(q)])) # symbol of t(Q), for some reason fft() applies t(Q) instead of Q
  
  m <- (K+delta-1)^2
  W0 <- matrix(rnorm(m),ncol=K+delta-1)
  w <- c(W0)
  
  w <- fft(q)*fft(w)
  w <- fft(w, inverse = T)/m
  w <- Re(w)
  W <- matrix(w,ncol=K+delta-1)
  W <- W[1:K,1:K]
  return(W)
}

create_B_stationary <- function(mask){
# creates covariance of the 2D WN process averaged using the `mask`,
# OUTPUT: a matrix of size delta x delta, which is the non-zero part of the symbol of B
  delta <- dim(mask)[1]
  B <- array(0,c(delta,delta))
  denom <- sum(mask^2)
  for(a in 1:delta){
    for(b in 1:delta){
      B[delta-a+1,delta-b+1] <- sum(mask[1:a,1:b]*mask[(delta-a+1):delta,(delta-b+1):delta])/denom
    }
  }
  return(B)
}

brownian_cov <- function(K){
# creates the covariance of a Brownian sheet on a K x K grid
  B <- array(0,c(K,K))
  for(i in 1:K){
    for(j in 1:K){
      B[i,j] <- min(i,j)
    }
  }
  return(B)
}

fourier_cov <- function(K){
  rank <- K
  V <- array(0,c(K,rank))
  x <- (1:K)/(K+1)
  for(r in 1:rank){
    if(r %% 2 == 1) V[,r] <- sin(2*pi*x*(r %/% 2 + 1)) else V[,r] <- cos(2*pi*x*(r %/% 2))
  }
  return(V %*% diag(2^(rank:1)) %*% t(V))
}

legendre_covariance <- function(K,r){
# creates C of rank with Legender polynomials
  X1_1D <- shifted_legendre(r, K, N)
  return(t(X1_1D)%*%X1_1D/N)
}

epanechnik <- function(x) 3/4*(1-x^2)*I(-1 <= x)*I(x <= 1)

create_data <- function(seed=517, buf=0, N=300, K=100, A1="legendre", A2="legendre", mask = NULL, tau=0.5, snr=NULL){
# INPUT: seed - to be set for reproducibility
#        N, K - number of observations and grid size
#         buf - averaging buffer, delta=2*buf+1, while delta=1 corresponds to corrupted diagonal
#      A1, A2 - `legendre` or `brownian` or specified covariances for the separable part
#        mask - by default all-ones delta x delta matrix, can be set to `epanechnik` or `triangular`
#         tau - in [0,1] parameter weighting together the banded part and the stationary part
# OUTPUT:   Y, X, W - generated data, Y = X + W
#         A1, A2, B - covariances used to generate `Y`, `B` is the symbol only (delta x delta matrix)
#                   - `tau` is already incorporated so cov(Y) = A1 x A2 + B
  delta <- 2*buf+1
  set.seed(seed)
  ranks <- c(7,7) # ranks for the legendre covariances
  # prepare covariance matrices for simulation
  if(is.null(A1)) A1 <- "legendre"
  if(is.null(A1)) A2 <- "legendre"
  if(A1=="legendre"){
    X1_1D <- shifted_legendre(ranks[1], K, N)
    A1 <- t(X1_1D)%*%X1_1D/N
  }else if(A1=="brownian") A1 <- brownian_cov(K)
  else if(A1=="fourier") A1 <- fourier_cov(K)
  if(A2=="legendre"){
    X2_1D <- shifted_legendre(ranks[2], K, N)
    A2 <- t(X2_1D)%*%X2_1D/N
  } else if(A2=="brownian") A2 <- brownian_cov(K)
  else if(A2=="fourier") A2 <- fourier_cov(K)
  
  if(is.null(mask)) mask <- matrix(1, delta, delta)
  else if(mask=="epanechnik") mask <- outer(epanechnik((-buf:buf)/(buf+1)),epanechnik((-buf:buf)/(buf+1)))
  else if(mask=="triangular") mask <- outer( (buf+1)-abs(-buf:buf) , (buf+1)-abs(-buf:buf))
  else{
    mask <- array(0,c(delta,delta))
    for(i in 1:delta){
      for(j in 1:delta){
        mask[i,j] <- (-1)^(i+j)
      }
    }
  }
  if(is.null(snr)){ # covariances are all standardized to have norm one
    A1 <- A1/sum(diag(A1))          
    A2 <- A2/sum(diag(A2))
    mask <- mask/frobenius(mask)/K
    B <- create_B_stationary(mask)/K^2
  }else{ # covariances standardized to have certain signal-to-noise ratio w.r.t. Frobenius norm
    A1 <- A1/frobenius(A1)           
    A2 <- A2/frobenius(A2)
    B <- create_B_stationary(mask)
    mask <- mask/frobenius(mask)/sqrt(snr*frobenius_stationary(K,B))
    B <- B/frobenius_stationary(K,B)/snr
  }
  
  A1_root <- rootmatrix(A1)
  A2_root <- rootmatrix(A2)

  W <- X <- array(0, dim = c(N, K, K))
  for (n in 1:N){
    X[n,,] <- t(A1_root)%*%matrix(rnorm(K^2), K)%*% A2_root
    W[n,,] <- create_banded_datum(mask,K)
  }
  X <- sqrt(1-tau)*X
  W <- sqrt(tau)*W
  Y <- X + W
  A1 <- sqrt(1-tau)*A1
  A2 <- sqrt(1-tau)*A2
  B <- tau*B
  return(list(Y=Y,X=X,W=W,A1=A1,A2=A2,B=B))
}

# create_data <- function(seed=517, buf=0, N=300, K=100, A1="legendre", A2="legendre", mask = NULL, tau=0.5){
#   # INPUT: seed - to be set for reproducibility
#   #        N, K - number of observations and grid size
#   #         buf - averaging buffer, delta=2*buf+1, while delta=1 corresponds to corrupted diagonal
#   #      A1, A2 - `legendre` or `brownian` or specified covariances for the separable part
#   #        mask - by default all-ones delta x delta matrix, can be set to `epanechnik` or `triangular`
#   #         tau - in [0,1] parameter weighting together the banded part and the stationary part
#   # OUTPUT:   Y, X, W - generated data, Y = X + W
#   #         A1, A2, B - covariances used to generate `Y`, `B` is the symbol only (delta x delta matrix)
#   #                   - `tau` is already incorporated so cov(Y) = A1 x A2 + B
#   delta <- 2*buf+1
#   set.seed(seed)
#   ranks <- c(7,7) # ranks for the legendre covariances
#   # prepare covariance matrices for simulation
#   if(is.null(A1)) A1 <- "legendre"
#   if(is.null(A1)) A2 <- "legendre"
#   if(A1=="legendre"){
#     X1_1D <- shifted_legendre(ranks[1], K, N)
#     A1 <- t(X1_1D)%*%X1_1D/N
#   }else if(A1=="brownian") A1 <- brownian_cov(K)
#   A1 <- A1/sum(diag(A1))           # covariances are all standardized to have trace one
#   A1_root <- rootmatrix(A1)
#   if(A2=="legendre"){
#     X2_1D <- shifted_legendre(ranks[2], K, N)
#     A2 <- t(X2_1D)%*%X2_1D/N
#   } else if(A2=="brownian") A2 <- brownian_cov(K)
#   A2 <- A2/sum(diag(A2))
#   A2_root <- rootmatrix(A2)
#   if(is.null(mask)) mask <- matrix(1, delta, delta)
#   else if(mask=="epanechnik") mask <- outer(epanechnik((-buf:buf)/(buf+1)),epanechnik((-buf:buf)/(buf+1)))
#   else if(mask=="triangular") mask <- outer( (buf+1)-abs(-buf:buf) , (buf+1)-abs(-buf:buf))
#   else{
#     mask <- array(0,c(delta,delta))
#     for(i in 1:delta){
#       for(j in 1:delta){
#         mask[i,j] <- (-1)^(i+j)
#       }
#     }
#   }
#   B <- create_B_stationary(mask)/K^2
#   
#   W <- X <- array(0, dim = c(N, K, K))
#   for (n in 1:N){
#     X[n,,] <- t(A1_root)%*%matrix(rnorm(K^2), K)%*% A2_root
#     W[n,,] <- create_banded_datum(mask,K)
#   }
#   X <- sqrt(1-tau)*X
#   W <- sqrt(tau)*W
#   Y <- X + W
#   A1 <- sqrt(1-tau)*A1
#   A2 <- sqrt(1-tau)*A2
#   B <- tau*B
#   return(list(Y=Y,X=X,W=W,A1=A1,A2=A2,B=B))
# }

################
### Plotting ### 
################

papcol <- function(i){
  bar <- c("forestgreen","goldenrod1","blue","blue","brown", "violet","red")
  return(bar[i])
}

pappch <- function(i){
  bar <- c(1,6,4,4,0,5,2)
  return(bar[i])
}

paplty <- function(i){
  bar <- c(4,2,3,3,1,5,6)
  return(bar[i])
}

