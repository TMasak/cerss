# require(lattice)
require(reshape2)
require(fdapace)
require(orthopolynom)
require(np)
# require(np)
# require(rgl)

estimate_sparse <- function(X,WN_flag=F,maxiter=2,B=NULL,bwidths=c(NULL,NULL),tol=1e-6,sparse=F){
# INPUT: X - data as a N x K1 x K2 array with NAs for missing values or a list of matrices with 3 columns (if sparse=T)
#  maxiter - no. of iterations to be computed
#        B - a candidate spatial covariance used as a starting point of the algorithm
#      tol - tolerance for the stopping criterion
#  bwidths - bandwidths for temporal and spatial smoothing. If not provided, CV is used in every step.
#  WN_flag - if true, a white noise errors are assumed and diagonals are discarded before every smoothign step
# OUTPUT: a list containing the estimated covariances A (temporal) and B (spatial)
# NOTE that data are assumed centered and only sparse=F is handled for now
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  if(length(B) == 0){
    B <- array(1, c(K2,K2))
  }
  A <- array(0,c(K1,K1))
  B_old <- array(0,c(K2,K2))
  eps <- 1
  iter <- 0
  while(iter < maxiter && eps > tol){
    A <- T1smooth(X,B,N,K1,K2,bwidths[1],WN_flag) # estimate temporal covariance given the spatial covariance
    A <- A$A/frobenius(A$A)                       # standardize to prevent mass escape
    B <- T2smooth(X,A,N,K1,K2,bwidths[2],WN_flag) # estimate spatial covariance given the temporal covariance
    B <- B$B
    iter <- iter +1
    eps <- frobenius(B-B_old)/frobenius(B_old)
    B_old <- B
    # print(iter)
  }
  # estimation of the error variance sigma^2
  if(WN_flag){
    Y <- X[,,]^2
    data <- melt(Y[1,,])
    data <- data[complete.cases(data),]
    for(n in 2:N){
      data <- rbind(data,melt(Y[n,,]))
      data <- data[complete.cases(data),]
    }
    Diag <- outer(diag(A),diag(B))
    obsGrid <- sort(unique(c(unlist(data[,1]))))
    pom <- as.matrix(data[,c(1,2)])
    attr(pom,"dimnames") <- NULL
    rcov <- list('tPairs' = pom, 'cyy' = data$value, 'cxxn' = data$value, 'error' = FALSE)
    if(K1==K2){ # PACE CV currently requires K1==K2 ... TODO: generalize to allow K1!=K2
      GCV <- GCVLwls2DV2(obsGrid = obsGrid, regGrid = 1:K1,rcov=rcov,t=data$Var1,kern='epan')
      bandwidth <- GCV$h
      # print(c("bandwidth sigma2",bandwidth))
    } else bandwidth <- min(K1,K2)/4
    Sigmasurf <- Lwls2D(bw=bandwidth, xin=as.matrix(data[,1:2]), yin=data$value, xout1=1:K1, xout2=1:K2)
    sigma2 <- mean(Sigmasurf[floor(K1/4):floor(K1-K1/4),floor(K2/4):floor(K2-K2/4)]) -
              mean(Diag[floor(K1/4):floor(K1-K1/4),floor(K2/4):floor(K2-K2/4)])
    return(list(A=A,B=B,sigma2=sigma2))
  } else return(list(A=A,B=B))
}

estimate_sparse_sim <- function(X,WN_flag=F,Cov1,Cov2,maxiter=2,B=NULL,bwidths=c(NULL,NULL),tol=1e-6,sparse=F){
# for separable covariance simulations only
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  if(length(B) == 0){
    B <- array(1, c(K2,K2))
  }
  A <- array(0,c(K1,K1))
  B_old <- array(0,c(K2,K2))
  eps <- 1
  iter <- 0
  err <- rep(0,maxiter)
  while(iter < maxiter && eps > tol){
    A <- T1smooth(X,B,N,K1,K2,bwidths[1],WN_flag) # estimate temporal covariance given the spatial covariance
    bw1 <- A$bw
    A <- A$A/frobenius(A$A)                       # standardize to prevent mass escape
    B <- T2smooth(X,A,N,K1,K2,bwidths[2],WN_flag) # estimate spatial covariance given the temporal covariance
    bw2 <- B$bw
    B <- B$B
    iter <- iter +1
    eps <- frobenius(B-B_old)/frobenius(B_old)
    B_old <- B
    err[iter] <- frobenius( outer(A,B) - outer(Cov1,Cov2) )
  }
  return(list(err=err,iter=iter,bw1=bw1,bw2=bw2))
}

estimate_sparse_gneiting <- function(X,WN_flag=F,C,maxiter=2,B=NULL,bwidths=c(NULL,NULL),tol=1e-6,sparse=F){
# for Gneiting's covariance simulations only
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  if(length(B) == 0){
    B <- array(1, c(K2,K2))
  }
  A <- array(0,c(K1,K1))
  B_old <- array(0,c(K2,K2))
  eps <- 1
  iter <- 0
  err <- rep(0,maxiter)
  while(iter < maxiter && eps > tol){
    A <- T1smooth(X,B,N,K1,K2,bwidths[1],WN_flag) # estimate temporal covariance given the spatial covariance
    bw1 <- A$bw
    A <- A$A/frobenius(A$A)                       # standardize to prevent mass escape
    B <- T2smooth(X,A,N,K1,K2,bwidths[2],WN_flag) # estimate spatial covariance given the temporal covariance
    bw2 <- B$bw
    B <- B$B
    iter <- iter +1
    eps <- frobenius(B-B_old)/frobenius(B_old)
    B_old <- B
    err[iter] <- frobenius( aperm(outer(A,B),c(1,3,2,4)) - C )
  }
  return(list(err=err,iter=iter,bw1=bw1,bw2=bw2))
}

T1smooth <- function(X,B,N,K1,K2,bandwidth,WN_flag=F){
# partial inner product w.r.t. to the first argument and subsequent smoothing, i.e. returning A (temporal)
  qq <- quantile(c(abs(B)),0.05) # 0.05 is quite arbitrary
  B[abs(B) < qq] <- 0 # discard where there is too little information
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Win <- Res <- array(0,c(N,K1,K1))
  if(WN_flag) diag(B) <- 0
  # prepare data points and their weights for temporal scatterplot smoothing
  for(n in 1:N){
    Ind <- !is.na(X[n,,])
    W <- Ind %*% B^2 %*% t(Ind)
    # W[W==0] <- NA
    W[abs(W) < sqrt(.Machine$double.eps)] <- NA
    X_akt <- X[n,,]
    X_akt[is.na(X_akt)] <- 0
    Res_akt <- X_akt %*% B %*% t(X_akt)
    Res_akt[Res_akt == 0] <- NA
    Res[n,,] <- Res_akt/W
    Win[n,,] <- W
  }
  # transform data to format required by PACE local linear smoothers
  Data <- melt(Res[1,,])
  Data <- Data[complete.cases(Data),]
  Weights <- melt(Win[1,,])
  Weights <- Weights[complete.cases(Weights),]
  for(n in 2:N){
    Data <- rbind(Data,melt(Res[n,,]))
    Data <- Data[complete.cases(Data),]
    Weights <- rbind(Weights,melt(Win[n,,]))
    Weights <- Weights[complete.cases(Weights),]
  }
  obsGrid <- sort(unique(c(unlist(Data[,1]))))
  pom <- as.matrix(Data[,c(1,2)])
  attr(pom,"dimnames") <- NULL
  rcov <- list('tPairs' = pom, 'cyy' = Data$value, 'cxxn' = Data$value, 'error' = FALSE)
  # use CV from the PACE package if bandwidth was not provided by user
  if(is.null(bandwidth)){
    GCV <- GCVLwls2DV2(obsGrid = obsGrid, regGrid = 1:K1,rcov=rcov,t=Data$Var1,kern='epan',Win = Weights$value)
    bandwidth <- GCV$h
    # print(c("Bandwidth A",bandwidth))
  }
  # use local linear smoother from the PACE package
  A <- Lwls2D(bw=bandwidth, xin=as.matrix(Data[,1:2]), yin=Data$value, xout1=1:K1, xout2=1:K1,win = Weights$value)
  return(list(A=A,bw=bandwidth))
}

T2smooth <- function(X,A,N,K1,K2,bandwidth,WN_flag=F){
# partial inner product w.r.t. to the second argument and subsequent smoothing, i.e. returning B (spatial)
  qq <- quantile(c(abs(A)),0.05) # 0.05 is quite arbitrary
  A[abs(A) < qq] <- 0 # discard where there is too little information
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Win <- Res <- array(0,c(N,K2,K2))
  if(WN_flag) diag(A) <- 0
  for(n in 1:N){
    Ind <- !is.na(X[n,,])
    W <- t(Ind) %*% A^2 %*% Ind
    # W[W==0] <- NA
    W[abs(W) < sqrt(.Machine$double.eps)] <- NA
    X_akt <- X[n,,]
    X_akt[is.na(X_akt)] <- 0
    Res_akt <- t(X_akt) %*% A %*% X_akt
    Res_akt[Res_akt == 0] <- NA
    Res[n,,] <- Res_akt/W
    Win[n,,] <- W
  }
  Data <- melt(Res[1,,])
  Data <- Data[complete.cases(Data),]
  Weights <- melt(Win[1,,])
  Weights <- Weights[complete.cases(Weights),]
  for(n in 2:N){
    Data <- rbind(Data,melt(Res[n,,]))
    Data <- Data[complete.cases(Data),]
    Weights <- rbind(Weights,melt(Win[n,,]))
    Weights <- Weights[complete.cases(Weights),]
  }
  obsGrid <- sort(unique(c(unlist(Data[,1]))))
  pom <- as.matrix(Data[,c(1,2)])
  attr(pom,"dimnames") <- NULL
  rcov <- list('tPairs' = pom, 'cyy' = Data$value, 'cxxn' = Data$value, 'error' = FALSE)
  if(is.null(bandwidth)){
    GCV <- GCVLwls2DV2(obsGrid = obsGrid, regGrid = 1:K2,rcov=rcov,t=Data$Var1,kern='epan',Win = Weights$value)
    bandwidth <- GCV$h
    # print(c("Bandwidth B",bandwidth))
  }
  B <- Lwls2D(bw=bandwidth, xin=as.matrix(Data[,1:2]), yin=Data$value, xout1=1:K2, xout2=1:K2,win = Weights$value)
  return(list(B=B,bw=bandwidth))
}

smooth4D <- function(X,bw1,bw2,WN_flag=T){
# estimates the covariance by local linear smoothing of a 4D scatterplot of raw covariances
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Times <- rep(0,2)
  if(!WN_flag){
    data <- melt(outer(X[1,,],X[1,,]))
    data <- data[complete.cases(data),]
    for(n in 2:N){
      data <- rbind(data,melt(outer(X[n,,],X[n,,])))
      data <- data[complete.cases(data),]
    }
    Xin <- as.matrix(data[,1:4])
    attr(Xin,"dimnames") <- NULL
    Xout <- melt(array(0,c(K1,K2,K1,K2)))
    Xout <- as.matrix(Xout[,1:4])
    attr(Xout,"dimnames") <- NULL
    # bw <- npregbw(xdat=Xin, ydat=data$value, exdat=Xout, regtype = "ll") # too expensive
    kre <- npreg(bws=c(bw1,bw2,bw1,bw2), xdat=Xin, ydat=data$value, exdat=Xout, regtype = "ll")
    C_hat <- array(kre$mean,c(K1,K2,K1,K2))
    return(list(C=C_hat))
  }else{
    t <- Sys.time() #
    Xx <- outer(X[1,,],X[1,,])
    Xx <- tensor2matrix(Xx)
    diag(Xx) <- NA                  # remove the raw covariances burdened by noise
    Xx <- matrix2tensor(Xx,K1,K2)
    data <- melt(Xx)
    data <- data[complete.cases(data),]
    for(n in 2:N){
      Xx <- outer(X[n,,],X[n,,])
      Xx <- tensor2matrix(Xx)
      diag(Xx) <- NA
      Xx <- matrix2tensor(Xx,K1,K2)
      data <- rbind(data,melt(Xx))
      data <- data[complete.cases(data),]
    }
    Times[1] <- difftime(Sys.time(),t,units="secs") #
    Xin <- as.matrix(data[,1:4])
    attr(Xin,"dimnames") <- NULL
    Xout <- melt(array(0,c(K1,K2,K1,K2)))
    Xout <- as.matrix(Xout[,1:4])
    attr(Xout,"dimnames") <- NULL
    t <- Sys.time() #
    kre <- npreg(bws=c(bw1,bw2,bw1,bw2), xdat=Xin, ydat=data$value, exdat=Xout, regtype = "ll")
    Times[2] <- difftime(Sys.time(),t,units="secs") #
    C_hat <- array(kre$mean,c(K1,K2,K1,K2))
    # noise level estimation
    Y <- X^2
    data <- melt(Y[1,,])
    data <- data[complete.cases(data),]
    for(n in 2:N){
      data <- rbind(data,melt(Y[n,,]))
      data <- data[complete.cases(data),]
    }
    Xin <- as.matrix(data[,1:2])
    attr(Xin,"dimnames") <- NULL
    Xout <- melt(array(0,c(K1,K2)))
    Xout <- as.matrix(Xout[,1:2])
    attr(Xout,"dimnames") <- NULL
    kre <- npreg(bws=c(bw1,bw2), xdat=Xin, ydat=data$value, exdat=Xout, regtype = "ll")
    Sigmasurf <- matrix(kre$mean,ncol=K2) - matrix(diag(tensor2matrix(C_hat)),ncol=K2)
    sigma2 <- mean(Sigmasurf[floor(K1/4):floor(K1-K1/4),floor(K2/4):floor(K2-K2/4)])
    return(list(C=C_hat,sigma2=sigma2,Times = Times))
  }
}

smooth4Dv1 <- function(X,bw1,bw2,WN_flag=T){
# estimates the covariance by local linear smoothing - only for timing purposes
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Times <- rep(0,2)
    t <- Sys.time() #
    Xx <- outer(X[1,,],X[1,,])
    Xx <- tensor2matrix(Xx)
    diag(Xx) <- NA                  # remove the raw covariances burdened by noise
    Xx <- matrix2tensor(Xx,K1,K2)
    data <- melt(Xx)
    data <- data[complete.cases(data),]
    for(n in 2:N){
      Xx <- outer(X[n,,],X[n,,])
      Xx <- tensor2matrix(Xx)
      diag(Xx) <- NA
      Xx <- matrix2tensor(Xx,K1,K2)
      data <- rbind(data,melt(Xx))
      data <- data[complete.cases(data),]
    }
    Times[1] <- difftime(Sys.time(),t,units="secs") #
    Xin <- as.matrix(data[,1:4])
    attr(Xin,"dimnames") <- NULL
    Xout <- melt(array(0,c(K1,K2,K1,K2)))
    Xout <- as.matrix(Xout[,1:4])
    attr(Xout,"dimnames") <- NULL
    t <- Sys.time() #
    kre <- npreg(bws=c(bw1,bw2,bw1,bw2), xdat=Xin, ydat=data$value, exdat=Xout, regtype = "ll")
    Times[2] <- difftime(Sys.time(),t,units="secs") #
    C_hat <- array(kre$mean,c(K1,K2,K1,K2))
    return(list(C=C_hat,Times = Times))
}

smooth4Dv2 <- function(X,WN_flag=T){
# estimates the covariance by loess - only for timing purposes
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Times <- rep(0,2)
    t <- Sys.time() #
    Xx <- outer(X[1,,],X[1,,])
    Xx <- tensor2matrix(Xx)
    diag(Xx) <- NA                  # remove the raw covariances burdened by noise
    Xx <- matrix2tensor(Xx,K1,K2)
    data <- melt(Xx)
    data <- data[complete.cases(data),]
    for(n in 2:N){
      Xx <- outer(X[n,,],X[n,,])
      Xx <- tensor2matrix(Xx)
      diag(Xx) <- NA
      Xx <- matrix2tensor(Xx,K1,K2)
      data <- rbind(data,melt(Xx))
      data <- data[complete.cases(data),]
    }
    Times[1] <- difftime(Sys.time(),t,units="secs") #
    Xin <- as.matrix(data[,1:4])
    attr(Xin,"dimnames") <- NULL
    Xout <- melt(array(0,c(K1,K2,K1,K2)))
    Xout <- as.matrix(Xout[,1:4])
    attr(Xout,"dimnames") <- NULL
    t <- Sys.time() #
    res <- loess(value~Var1+Var2+Var3+Var4, span=0.25, data=data, degree=1)
    C_hat_vals <- predict(res, newdata=Xout)
    C_hat <- array(C_hat_vals, c(K1,K2,K1,K2))
    Times[2] <- difftime(Sys.time(),t,units="secs") #
    return(list(C=C_hat,Times = Times))
}

### non-pooled observations, squared weights

estimate_sparse_nonpooled <- function(X,WN_flag=F,maxiter=2,B=NULL,bwidths=c(NULL,NULL),tol=1e-6,sparse=F){
  # INPUT: X - data as a N x K1 x K2 array with NAs for missing values or a list of matrices with 3 columns (if sparse=T)
  #  maxiter - no. of iterations to be computed
  #        B - a candidate spatial covariance used as a starting point of the algorithm
  #      tol - tolerance for the stopping criterion
  #  bwidths - bandwidths for temporal and spatial smoothing. If not provided, CV is used in every step.
  #  WN_flag - if true, a white noise errors are assumed and diagonals are discarded before every smoothign step
  # OUTPUT: a list containing the estimated covariances A (temporal) and B (spatial)
  # NOTE that data are assumed centered and only sparse=F is handled for now
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  if(length(B) == 0){
    B <- array(1, c(K2,K2))
  }
  A <- array(0,c(K1,K1))
  B_old <- array(0,c(K2,K2))
  eps <- 1
  iter <- 0
  while(iter < maxiter && eps > tol){
    A <- T1smooth4(X,B,N,K1,K2,bwidths[1],WN_flag) # estimate temporal covariance given the spatial covariance
    A <- A/frobenius(A)                           # standardize to prevent mass escape
    B <- T2smooth4(X,A,N,K1,K2,bwidths[2],WN_flag) # estimate spatial covariance given the temporal covariance
    iter <- iter +1
    eps <- frobenius(B-B_old)/frobenius(B_old)
    B_old <- B
  }
  # estimation of the error variance sigma^2
  # if(WN_flag){
  #   Y <- X[,,]^2
  #   data <- melt(Y[1,,])
  #   data <- data[complete.cases(data),]
  #   for(n in 2:N){
  #     data <- rbind(data,melt(Y[n,,]))
  #     data <- data[complete.cases(data),]
  #   }
  #   Diag <- outer(diag(A),diag(B))
  #   obsGrid <- sort(unique(c(unlist(data[,1]))))
  #   pom <- as.matrix(data[,c(1,2)])
  #   attr(pom,"dimnames") <- NULL
  #   rcov <- list('tPairs' = pom, 'cyy' = data$value, 'cxxn' = data$value, 'error' = FALSE)
  #   if(K1==K2){ # PACE CV currently requires K1==K2 ... TODO: generalize to allow K1!=K2
  #     GCV <- GCVLwls2DV2(obsGrid = obsGrid, regGrid = 1:K1,rcov=rcov,t=data$Var1,kern='epan')
  #     bandwidth <- GCV$h
  #     # print(c("bandwidth sigma2",bandwidth))
  #   } else bandwidth <- min(K1,K2)/4
  #   Sigmasurf <- Lwls2D(bw=bandwidth, xin=as.matrix(data[,1:2]), yin=data$value, xout1=1:K1, xout2=1:K2)
  #   sigma2 <- mean(Sigmasurf[floor(K1/4):floor(K1-K1/4),floor(K2/4):floor(K2-K2/4)]) -
  #     mean(Diag[floor(K1/4):floor(K1-K1/4),floor(K2/4):floor(K2-K2/4)])
  #   return(list(A=A,B=B,sigma2=sigma2))
  # } else # return below here
  return(list(A=A,B=B))
}

T1smooth4 <- function(X,B,N,K1,K2,bandwidth,WN_flag=F){
  # partial inner product w.r.t. to the first argument and subsequent smoothing, i.e. returning A (temporal)
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  if(WN_flag) diag(B) <- 0
  # prepare data points and their weights for temporal scatterplot smoothing
  Data <- melt(outer(X[1,,],X[1,,]))
  Data <- Data[complete.cases(Data),]
  Weights <- Data
  for(m in 1:nrow(Data)){
    Data$value[m] <- Data$value[m]/B[Data$Var2[m],Data$Var4[m]]
    Weights$value[m] <- B[Data$Var2[m],Data$Var4[m]]^2
  }
  for(n in 2:N){
    Y <- melt(outer(X[n,,],X[n,,]))
    Y <- Y[complete.cases(Y),]
    W <- Y
    for(m in 1:nrow(Y)){
      Y$value[m] <- Y$value[m]/B[Y$Var2[m],Y$Var4[m]]
      W$value[m] <- B[Y$Var2[m],Y$Var4[m]]^2
    }
    Data <- rbind(Data,Y)
    Weights <- rbind(Weights,W)
  }
  Data <- Data[,-c(2,4)]
  Weights <- Weights[,-c(2,4)]
  # # pooling observations before smoothing?
  # pooled <- pool_obs(Data,Weights,K1)
  # Data <- pooled$Y
  # Weights <- pooled$W
  # #
  obsGrid <- sort(unique(c(unlist(Data[,1]))))
  pom <- as.matrix(Data[,c(1,2)])
  attr(pom,"dimnames") <- NULL
  rcov <- list('tPairs' = pom, 'cyy' = Data$value, 'cxxn' = Data$value, 'error' = FALSE)
  # use CV from the PACE package if bandwidth was not provided by user
  if(is.null(bandwidth)){
    GCV <- GCVLwls2DV2(obsGrid = obsGrid, regGrid = 1:K1,rcov=rcov,t=Data$Var1,kern='epan',Win = Weights$value)
    bandwidth <- GCV$h
    # print(c("Bandwidth A",bandwidth))
  }
  # use local linear smoother from the PACE package
  A <- Lwls2D(bw=bandwidth, xin=as.matrix(Data[,1:2]), yin=Data$value, xout1=1:K1, xout2=1:K1,win = Weights$value)
  return(A)
}

T2smooth4 <- function(X,A,N,K1,K2,bandwidth,WN_flag=F){
  # partial inner product w.r.t. to the second argument and subsequent smoothing, i.e. returning B (spatial)
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  if(WN_flag) diag(A) <- 0
  # prepare data points and their weights for temporal scatterplot smoothing
  Data <- melt(outer(X[1,,],X[1,,]))
  Data <- Data[complete.cases(Data),]
  Weights <- Data
  for(m in 1:nrow(Data)){
    Data$value[m] <- Data$value[m]/A[Data$Var1[m],Data$Var3[m]]
    Weights$value[m] <- A[Data$Var1[m],Data$Var3[m]]^2
  }
  for(n in 2:N){
    Y <- melt(outer(X[n,,],X[n,,]))
    Y <- Y[complete.cases(Y),]
    W <- Y
    for(m in 1:nrow(Y)){
      Y$value[m] <- Y$value[m]/A[Y$Var1[m],Y$Var3[m]]
      W$value[m] <- A[Y$Var1[m],Y$Var3[m]]^2
    }
    Data <- rbind(Data,Y)
    Weights <- rbind(Weights,W)
  }
  Data <- Data[,-c(1,3)]
  Weights <- Weights[,-c(1,3)]
  # # pooling observations before smoothing?
  # pooled <- pool_obs(Data,Weights,K2)
  # Data <- pooled$Y
  # Weights <- pooled$W
  # #
  obsGrid <- sort(unique(c(unlist(Data[,1]))))
  pom <- as.matrix(Data[,c(1,2)])
  attr(pom,"dimnames") <- NULL
  rcov <- list('tPairs' = pom, 'cyy' = Data$value, 'cxxn' = Data$value, 'error' = FALSE)
  if(is.null(bandwidth)){
    GCV <- GCVLwls2DV2(obsGrid = obsGrid, regGrid = 1:K2,rcov=rcov,t=Data$Var2,kern='epan',Win = Weights$value)
    bandwidth <- GCV$h
    # print(c("Bandwidth B",bandwidth))
  }
  B <- Lwls2D(bw=bandwidth, xin=as.matrix(Data[,1:2]), yin=Data$value, xout1=1:K2, xout2=1:K2,win = Weights$value)
  return(B)
}

###############
### fdapace ### some functions from the PACE package that needed to be modified a bit
###############

GCVLwls2DV2 <- function(obsGrid, regGrid, ngrid=NULL, dataType=rcov$dataType, error=rcov$error, kern, 
                        rcov, h0=NULL, verbose=FALSE, CV=FALSE, t, useBW1SE = FALSE, Win=NULL) {
  
  # TODO:? get the residual values only within truncated regGrid
  
  # Returns: a list of length 2, containing the optimal bandwidth and the gcv score.
  # obsGrid: observation points. 
  # ngrid: I think this should not be used in the gcv function.
  # CV: whether to use CV rather than GCV. Default to FALSE as not using CV. If CV is used use an integer value to specify the number of cross-validation folds. 
  
  # This function computes the optimal bandwidth choice for the covariance surface. 
  # based on the  GCV method by pooling the longitudinal data together. 
  # verbose is unused for now
  # this follows exactly the matlab 2D gcv selector.
  
  ## Get minimal bandwidth and range
  r <- diff(range(obsGrid)) * sqrt(2) # sqrt(2) because the window is circular.
  minBW <- fdapace:::GetMinb(t, dataType=rcov$dataType, obsGrid=obsGrid)
  
  if (missing(h0)) {
    h0 <- minBW
  }
  
  if (is.null(h0)){
    stop('the data is too sparse, no suitable bandwidth can be found! Try Gaussian Kernel instead!\n')
  }
  
  if (kern == 'gauss') {
    h0 = h0 * 0.2;
  }
  
  ## Get Candidate Bandwidths
  h0 <- min(h0, r/4)
  if (h0 < r/4) {    
    q <- (r / (4 * h0)) ^ (1/9)
  } else if (h0 < r/2) {
    q <- (r / (2 * h0)) ^ (1/9)
  } else if (h0 < r) {
    q <- (r / h0) ^ (1/9)
  } else {
    stop('Data is too sparse. The minimal bandwidth is the range of data')
  }
  bw <- (q ^ seq(0,9,length.out = 10)) * h0 # from h0 to r / 4
  
  
  ## Set GCV/CV options
  opth <- h0
  
  leave <- FALSE
  iter <- 0
  maxIter <- 1
  if (CV != FALSE) {
    # We partition the raw covariance rather than partition the individuals.
    fold <- CV
    partition <- CreateFolds(1:nrow(rcov$tPairs), k=fold)
  }
  
  minBWInvalid <- FALSE
  while (!leave && iter < maxIter) {
    if (minBWInvalid){
      minBW <- bw[1]
    } 
    
    #Scores <- rep(Inf, length(bw))
    Scores <- matrix(Inf, nrow = length(bw), ncol = 2); colnames(Scores) <- c('SUM','SD');
    # try the bandwidths large to small in order to save time due to sparseness in the windows.
    for (i in rev(seq_along(bw))) {
      h <- bw[i]
      
      if (class(rcov) == 'BinnedRawCov') {
        if (CV == FALSE) # GCV
          Scores[i,'SUM'] <- getGCVscoresV2(h, kern, rcov$tPairs, rcov$meanVals, win=rcov$count, regGrid=regGrid, RSS=rcov$RSS, verbose=verbose)
        else # CV
          Scores[i,] <- getCVscoresV2(partition, h, kern, rcov$tPairs, rcov$meanVals, win=rcov$count, regGrid=regGrid, RSS=rcov$RSS, verbose=verbose)
      } else {
        if (CV == FALSE) # GCV
          Scores[i,'SUM'] <- getGCVscoresV2(h, kern, rcov$tPairs, rcov$cxxn, regGrid=regGrid, verbose=verbose)
        else # CV
          Scores[i,] <- getCVscoresV2(partition, h, kern, rcov$tPairs, rcov$cxxn,win = Win, regGrid=regGrid, verbose=verbose)
      }
      
      if (is.infinite(Scores[i,'SUM'])) { 
        minBWInvalid <- TRUE
        if (i < length(bw)) {
          if (minBWInvalid) {
            minBW <- bw[i + 1]
            minBWInvalid <- FALSE
          }
        }
        break; # This will help break out of the loop if the BW is too small to make sense
      }
    } 
    
    #optInd <- which.min(Scores)
    #opth <- bw[optInd]
    #optgcv <- Scores[optInd]
    
    if(is.infinite(min(Scores[,'SUM']))){
      opth <- max(bw)
      optgcv <- Inf
      # } else if( sum(!is.infinite(Scores)) >= 2 ){ # Given we have at least two points we can fit "something"
      # nonInf = which(!is.infinite(Scores));
      # costSpline = spline(bw[nonInf], Scores[nonInf])
      # opth = costSpline$x[which.min(costSpline$y)]
      # optgcv = min(costSpline$y)
    } else {
      if(useBW1SE){
        ind = max(which( (Scores[,'SUM']/fold)  < (min(Scores[,'SUM'])/fold) + Scores[which.min(Scores[,'SUM']),'SD']/sqrt(fold) ))
        opth <- bw[ind]
        optgcv <- Scores[ind,'SUM']
      } else {
        ind <- which.min(Scores[,'SUM'])
        opth <- bw[ind]
        optgcv <- Scores[ind,'SUM']
      }
    } 
    
    ## Check that what we found is coherent.
    if (opth >= r - 1e-12) {
      minBW <- r
      leave <- TRUE
      stop('Data is too sparse. The optimal bandwidth equals to the range of input time points. Try Gaussian kernel.')
    }
    if ( (abs(opth - max(bw)) > 1e-12) && !is.infinite(optgcv))
      leave <- TRUE            
    else if (is.infinite(optgcv)) {
      if (verbose)
        warning('Data is too sparse, retry with larger bandwidths!')
      h0 <- max(bw) * 1.01
    } else if ( (abs(opth - max(bw)) > 1e-12) ) {
      warning('Optimal bandwidth not found in the candidate bandwidths. Retry with larger bandwidths')
      h0 <- max(bw)
    }
    if (!leave) {
      newr <- seq(0.5, 1, by=0.05) * r # ??? this can be quite slow
      ind <- which(newr > h0)[1]
      q <- (newr[ind] / h0) ^ (1/9)
      bw <- q ^ (0:9) * h0
      if (verbose) {
        message('New bwuserCov candidates:\n')
        print(bw)
      }
      iter <- iter + 1
    }
  } # The "while (!leave && iter < maxIter) ..." end here
  
  ret <- list(h=opth, gcv=optgcv, minBW=minBW)
  if (CV != FALSE)
    names(ret)[2] <- 'cv'
  
  return(ret)
  
}


getGCVscoresV2 <- function(bw, kern, xin, yin, win=NULL, regGrid, RSS=NULL, verbose=FALSE) {
  # ...: passed on to Lwls2D
  # RSS: for implementing GCV of binned rcov.
  # browser() 
  if (is.null(win))
    win <- rep(1, length(yin))
  
  fit <- tryCatch(Lwls2D(bw, kern, xin=xin, yin=yin, win=win, xout1=regGrid, xout2=regGrid), error=function(err) {
    if (verbose) {
      message('Invalid bandwidth. Try enlarging the window size.\n')
    }
    return(Inf)
  })
  
  # Catch
  if (is.infinite(fit[1]))
    return(Inf)
  
  # workaround for degenerate case.
  if (any(is.nan(fit)))
    return(Inf)
  
  obsFit <- fdapace:::interp2lin(as.numeric(regGrid), as.numeric(regGrid), fit, as.numeric(xin[, 1]), as.numeric(xin[, 2]))
  
  # residual sum of squares
  res <- sum((yin - obsFit) ^ 2 * win)
  if (!is.null(RSS))
    res <- res + sum(RSS)
  
  # kernel at x=0
  k0 <- KernelAt0(kern)
  N <- sum(win)
  r <- diff(range(xin[, 1]))
  bottom <- max(1 - 3 * (1 / N) * (r * k0 / bw)^2, 0)
  GCV <- res / bottom^2
  
  return(GCV)
}

# k-fold CV
# partition: a list of testset observation indices, returned by caret::createFolds
# ...: passed on to Lwls2D
getCVscoresV2 <- function(partition, bw, kern, xin, yin, win=NULL, regGrid, RSS=NULL, verbose=FALSE) {
  
  if (is.null(win))
    win <- rep(1, length(yin))
  
  n <- length(yin)
  
  # browser()
  cvSubSum <- sapply(partition, function(testSet) {
    # browser()
    fit <- tryCatch(Lwls2D(bw, kern, xin=xin, yin=yin, win=win, xout1=regGrid, xout2=regGrid, subset=-testSet), error=function(err) {
      if (verbose) {
        message('Invalid bandwidth. Try enlarging the window size.\n')
      }
      return(Inf)
    })
    
    # Catch
    if (is.infinite(fit[1]))
      return(Inf)
    
    # workaround for degenerate case.
    if (any(is.nan(fit)))
      return(Inf)
    
    obsPred <- interp2lin(regGrid, regGrid, fit, xin[testSet, 1], xin[testSet, 2])
    tmpSum <- sum((yin[testSet] - obsPred) ^ 2 * win[testSet])
    
    # residual sum of squares
    if (!is.null(RSS))
      tmpSum <- tmpSum + sum(RSS[testSet])
    
    return(tmpSum)
  })
  
  return(c(sum(cvSubSum), sd(cvSubSum)))
}

KernelAt0 <- function(kern) {
  if (kern == 'quar')
    k0 <- 0.9375
  else if (kern == 'epan')
    k0 <- 0.75
  else if (kern == 'rect')
    k0 <- 0.5
  else if (kern == 'gausvar')
    k0 <- 0.498677850501791
  else if (kern == 'gauss')
    k0 <- 0.398942280401433
  else
    stop('Unknown kernel')
  
  return(k0)
}

#################
### utilities ###
#################

frobenius <- function(X){
  return(sqrt(sum(X^2)))
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

rootmatrix <- function(x){
  # Calculates "square-root" of a given matrix, taken from the elasticnet package. 
  x.eigen <- eigen(x)
  d <- x.eigen$values
  d <- (d + abs(d))/2
  v <- x.eigen$vectors
  return(v %*% diag(sqrt(d)) %*% t(v))
}

plotni3D <- function(Data, A=NULL, Rainbow=rainbow(20), psize=2){
# plots scatterplot smoothing in 3D, z-axis has no meaning here, values rescaled for plotting purposes
# Input: Data is the same data.frame that is passed to Lwls2D() of 'fdapace' and A is the output of Lwls2s()
#        Rainbow - vector of colors for heat-colored surface plot. E.g.: "green" or rainbow(20)
  # first create the 3D plot and plot the scattered points
  K1 <- max(Data$Var1)
  K2 <- max(Data$Var2)
  skalar <- max(K1,K2)/(max(Data$value,na.rm=T)-min(Data$value,na.rm=T))/1.5 # NA because of NA's
  rgl.open()# Open a new RGL device
  rgl.bg(color = "white") # Setup the background color
  rgl.points(Data$Var1,Data$value*skalar,Data$Var2, color = "blue", size = psize) # Scatter plot
  rgl.bbox(color=c("gray100","black"),xat=c(1,seq(5,K1,by=5)), zat = c(1,seq(5,K2,by=5))) 
  rgl.viewpoint( theta = -120, phi = 20)
  # now plot the fitted surface
  if(!is.null(A)){
    K1 <- dim(A)[1]
    K2 <- dim(A)[2] # K1 and K2 are mostly the same, but in theory A can be fitted by Lwls2D() on non-square grid
    Lines <- array(0,c(2*K1*(K2-1),3)) # x,y,z coordinates, x \in 1:K1, y in 1:K2, z is surface value
    Colors <- rep("",2*K1*(K2-1))
    colorscore <- 0
    Range <- diff(range(A))
    for(i in 1:K1){
      for(j in 1:(K2-1)){
        Lines[(i-1)*2*(K2-1)+2*j-1,1] <- i
        Lines[(i-1)*2*(K2-1)+2*j-1,2] <- j
        Lines[(i-1)*2*(K2-1)+2*j-1,3] <- A[i,j]
        Lines[(i-1)*2*(K2-1)+2*j,1] <- i
        Lines[(i-1)*2*(K2-1)+2*j,2] <- j+1
        Lines[(i-1)*2*(K2-1)+2*j,3] <- A[i,j+1]
        colorscore <- floor(((A[i,j] + A[i,j+1])/2-min(A))/(Range+1e-8)*length(Rainbow))+1
        Colors[(i-1)*2*(K2-1)+2*j-1] <- Rainbow[colorscore]
        Colors[(i-1)*2*(K2-1)+2*j-1+1] <- Rainbow[colorscore]
      }
    }
    rgl.lines(Lines[,1],Lines[,3]*skalar,Lines[,2], color = Colors)
    Lines <- array(0,c(2*K2*(K1-1),3)) # x,y,z coordinates, x \in 1:K1, y in 1:K2, z is surface value
    Colors <- rep("",2*K1*(K2-1))
    colorscore <- 0
    Range <- diff(range(A))
    for(j in 1:K2){
      for(i in 1:(K1-1)){
        Lines[(j-1)*2*(K1-1)+2*i-1,1] <- i
        Lines[(j-1)*2*(K1-1)+2*i-1,2] <- j
        Lines[(j-1)*2*(K1-1)+2*i-1,3] <- A[i,j]
        Lines[(j-1)*2*(K1-1)+2*i,1] <- i+1
        Lines[(j-1)*2*(K1-1)+2*i,2] <- j
        Lines[(j-1)*2*(K1-1)+2*i,3] <- A[i+1,j]
        colorscore <- floor(((A[i+1,j] + A[i,j])/2-min(A))/(Range+1e-8)*length(Rainbow))+1
        Colors[(j-1)*2*(K1-1)+2*i-1] <- Rainbow[colorscore]
        Colors[(j-1)*2*(K1-1)+2*i-1+1] <- Rainbow[colorscore]
      }
    }
    rgl.lines(Lines[,1],Lines[,3]*skalar,Lines[,2], color = Colors)
  }
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
  lb <- 1
  ub <- 2^(rank-1)
  lambda <- sort(c(3,seq(lb,ub,length.out=rank)),decreasing = TRUE)
  C <- C*lambda
  pol_val <- t(pol_val%*%C)
  pol_val <- t(t(pol_val)-apply(pol_val,2,mean))
  return(pol_val)
}

get_cov <- function(K, covtype="brownian",N=NULL, rank=NULL){
# K is the grid size, returns a covariance matrix of size K x K
# N and rank are for the empirical legendre covariance
  if(covtype=="brownian"){
    C <- array(0,c(K,K))
    for(i in 1:K){
      for(j in 1:K){
        C[i,j] <- min(i,j)
      }
    }
  } else if(covtype=="legendre"){
    if(is.null(N)) N <- 100
    if(is.null(rank)) rank <- 4
    X1_1D <- shifted_legendre(rank, K, N)
    C <- t(X1_1D)%*%X1_1D/N
  } else if(covtype=="fourier"){
    if(is.null(rank)) rank <- K
    V <- array(0,c(K,rank))
    x <- (1:K)/(K+1)
    for(r in 1:rank){
      if(r %% 2 == 1) V[,r] <- sin(2*pi*x*(r %/% 2 + 1)) else V[,r] <- cos(2*pi*x*(r %/% 2))
    }
    C <- V %*% diag(2^(rank:1)) %*% t(V)
  }
  return(C)
}

simulate_data <- function(N, perc, A, B, vary=F, sparse=F, WN_flag=F, save_as=NULL){
# INPUT: N, K1, K2 as usual - no. of surfaces, grid size in temporal dimension, grid size in spatial dimension
#        perc - percentage of observations per a single surface (i.e. K1*K2*perc/100 points per surface)
#         A,B - temporal and spatial covariances
#        vary - if False, the same no. of points is observed for every surface (depending on perc)
#             - if True, the no. of points per surface is uniform from 1 to a value given by 2*perc
#      sparse - how data is stored. If False, dense pattern is returned, i.e. a N x K1 x K2 array with (a lot of ) NAs.
#               If True, a sparse patter is returned, i.e. a list of N entries, every entry being a matrix of 3 columns
#               for temporal axis value, spatial axis value, and the value of the surface
#     WN_flag - whether a White Noise errors should be added
#     save_as - if a string such as "data100" is provided, the data are saved in the sparse pattern into "data100.txt"
# Output: data set Y in the format specified by 'sparse'
  K1 <- dim(A)[1]
  K2 <- dim(B)[2]
  # simulate data
  A_root <- rootmatrix(A)
  B_root <- rootmatrix(B)
  X <- array(0, dim = c(N, K1, K2))
  for (n in 1:N){
    X[n,,] <- t(A_root)%*%matrix(rnorm(K1*K2), K1)%*% B_root
  }
  if(WN_flag) X <- X + 0.05*array(rnorm(N*K1*K2),c(N,K1,K2))
  # subsample data
  for(n in 1:N){
    if(vary) n_points <- sample(1:(2*floor(K1*K2*perc/100)),1)
    else n_points <- floor(K1*K2*perc/100)
    inds <- sample(1:(K1*K2),n_points)
    # inds <- sample(c(1:(K1*K2),1:(K1*K2/2)),n_points) # systematic oversampling half of the domain
    X_akt <- X[n,,]
    X_akt[-inds] <- NA
    X[n,,] <- X_akt
  } 
  # if asked, transform data into the sparse format
  Y <- X
  if(sparse){
    X <- list()
    for(n in 1:N){
      X_akt <- melt(Y[n,,])
      X[[n]] <- X_akt[complete.cases(X_akt),]
    }
  }
  # if asked, save the data as a text file
  if(sparse && !is.null(save_as)){
    print("TODO")
  }
  return(X)
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
  # gamma <- 1/2  # these values here are fitted in the Genton paper.
  # alpha <- 0.834 
  # sigma2 <- 1 - 0.0415 
  # a <- 0.972 # 2500
  # const <- 0.00128 # 2500 # denoted 'c' in Bagchi & Dette (2017)
  # tau <- 1
  gamma <- 1 # these values were used in our simulation study
  alpha <- 1
  sigma2 <- 1
  a <- 1 # 2500
  const <- 1 # 2500 # denoted 'c' in Bagchi & Dette (2017)
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

run_sim <- function(N,K1,K2,perc,maxiter=20,seed=517,covtype="brownian",WN_flag=F,vary=F,sparse=F){
  # simulate data
  set.seed(seed)
  A <- get_cov(K1,covtype)
  A <- A/frobenius(A)
  B <- get_cov(K2,covtype)
  B <- B/frobenius(B)
  X <- simulate_data(N,perc,A,B,vary,sparse,WN_flag)
# T1smooth(X,B,N,K1,K2,F)
  Res <- estimate_sparse(X,maxiter,B=B,WN_flag=WN_flag,sparse=F)
  err <- frobenius( outer(Res$A,Res$B) - outer(A,B))
  return(err)
}

