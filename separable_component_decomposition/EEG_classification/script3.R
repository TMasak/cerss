# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("library_classification.R")

classify <- function(foldno){
  # load data
  alcoholic <- as.logical(c(rep(1,65),rep(0,44),rep(1,12),0))
  alcoholic <- alcoholic[-47]
  load("aData3.RData")
  # choose 5 surfaces for prediction
  set.seed(517)
  perm <- sample(1:121)
  Folds <- matrix(perm[1:120],ncol=5)
  Dat_test <- aData3[Folds[foldno,],,]
  Dat_train <- aData3[-Folds[foldno,],,]
  alc_train <- alcoholic[-Folds[foldno,]]
  alc_test <- alcoholic[Folds[foldno,]]
  # fit the covariance models
  mu0 <- apply(Dat_train[!alc_train,,],c(2,3),mean)
  mu1 <- apply(Dat_train[alc_train,,],c(2,3),mean)
  Dat_train[!alc_train,,] <- sweep(Dat_train[!alc_train,,], c(2,3), mu0)
  Dat_train[alc_train,,] <- sweep(Dat_train[alc_train,,], c(2,3), mu1)

  Res2 <- scd(Dat_train,R=2)
  Res1 <- Res2
  Res1$A <- Res1$A[1,,,drop=F]
  Res1$B <- Res1$B[1,,,drop=F]
  Res1$sigma <- Res1$sigma[1]
  
  n_reg <- 10 # how many different regularizers are we going to try out
  Res <- array(0,c(2,n_reg))
  Eps <- 1:n_reg/(n_reg+1)
  Eps <- sort(log(Eps)^2/2)
  
  for(j in 1:n_reg){
    print(j)
    psi1 <- ridge_psi(mu0,mu1,Res1,Eps[j])
    psi2 <- ridge_psi(mu0,mu1,Res2,Eps[j])
    for(k in 1:5){
      class <- LDA(mu0,mu1,Res1,Dat_test[k,,],psi1,F)
      Res[1,j] <- Res[1,j] + abs(alc_test[k]-class)
      class <- LDA(mu0,mu1,Res2,Dat_test[k,,],psi2,F)
      Res[2,j] <- Res[2,j] + abs(alc_test[k]-class)
    }
  }
  return(Res)
}

library(parallel)
Res <- mclapply(1:24,classify,mc.cores=25)
save(Res, file="res3ridge.RData")



