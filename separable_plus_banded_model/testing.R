setwd("C:/Users/Tomas/Documents/Skola/EPFL/Tex/PT_paper/JASA_resubmission/codes")
source("library.R")
library(ggplot2)
library(gridExtra)
library(reshape2)
library(scatterplot3d)
library(locfit)
library(covsep)

###################
### death rates ###
###################

load("mortality_data.RData")

testres <- empirical_bootstrap_test(data,2,4,"no") # 2,4=0.06; 4,2=0.07; 2,2=0.66; 4,4=0.076
testres2 <- empirical_bootstrap_test2(data,2,4,"no") # 2,4=0.367 (up to 0.42); 4,2=0.397; 1,1=

sdata <- data
for(n in 1:32){
  country <- data[n,,]
  country <- melt(country)
  smoo <- smooth.lf(x=as.matrix(country[,1:2]), y=country$value, kern = "epan",kt="prod",deg=1,alpha=c(0,2),maxk=2000)
  sdata[n,,] <- matrix(smoo$y, ncol=40)
}
testress <- empirical_bootstrap_test(sdata,2,4,"no") # 2,4=0.29;




empirical_bootstrap_test2 <- function (Data, L1 = 1, L2 = 1, studentize = "full", B = 1000,
          verbose = TRUE)
{
  N <- dim(Data)[1]
  d1 <- dim(Data)[2]
  d2 <- dim(Data)[3]
  # marginal.cov <- marginal_covariances(Data)
  proj.diff <- projected_differences2(Data, max(L1), max(L2),
                                     with.asymptotic.variances = TRUE)
  stat = array(NA, length(L1))
  for (k in 1:length(L1)) {
    l1 = L1[k]
    l2 = L2[k]
    if (studentize == "no") {
      stat[k] = sum(proj.diff$T.N[1:l1, 1:l2]^2)
    }
    else if (studentize == "diag" || studentize ==
             "full") {
      stat[k] = sum(renormalize_mtnorm(proj.diff$T.N[1:l1,
                                                     1:l2, drop = FALSE], proj.diff$sigma.left[1:l1,
                                                                                               1:l1, drop = FALSE], proj.diff$sigma.right[1:l2,
                                                                                                                                          1:l2, drop = FALSE], type = studentize)^2)
    }
    else {
      stop(paste0("unknown value for studentize: ",
                  studentize))
    }
  }
  if (verbose)
    cat("*** Empirical Bootstrap ***\n")
  boot.stat <- matrix(NA, B, length(L1))
  for (b in 1:B) {
    if (verbose) {
      if ((b%%floor(1 + B/200)) == 0) {
        cat("EB", round(100 * b/B, 1), "%--",
            sep = "")
      }
    }
    Data.boot <- Data[sample.int(N, N, replace = TRUE), ,
    ]
    proj.diff.boot = projected_differences2(Data.boot, max(L1),
                                           max(L2), with.asymptotic.variances = TRUE)
    for (k in 1:length(L1)) {
      l1 = L1[k]
      l2 = L2[k]
      if (studentize == "no") {
        boot.stat[b, k] <- sum((proj.diff.boot$T.N[1:l1,
                                                   1:l2] - proj.diff$T.N[1:l1, 1:l2])^2)
      }
      else if (studentize == "diag" || studentize ==
               "full") {
        boot.stat[b, k] = sum(renormalize_mtnorm(proj.diff.boot$T.N[1:l1,
                                                                    1:l2, drop = FALSE] - proj.diff$T.N[1:l1, 1:l2,
                                                                                                        drop = FALSE], proj.diff.boot$sigma.left[1:l1,
                                                                                                                                                 1:l1, drop = FALSE], proj.diff.boot$sigma.right[1:l2,
                                                                                                                                                                                                 1:l2, drop = FALSE], type = studentize)^2)
      }
      else {
        stop(paste0("unknown value for studentize: ",
                    studentize))
      }
    }
    rm(proj.diff.boot)
  }
  pvalues <- rep(NA, dim(boot.stat)[2])
  for (k in 1:length(L1)) {
    pvalues[k] <- mean(boot.stat[, k] > stat[k])
  }
  return(pvalues)
}

projected_differences2 <- function (Data, l1 = 1, l2 = 1, with.asymptotic.variances = TRUE)
{
  N <- dim(Data)[1]
  d1 <- dim(Data)[2]
  d2 <- dim(Data)[3]
  Data = sweep(Data, c(2, 3), apply(Data, c(2, 3), mean))
  marginal.covariances = estimate_separable(Data,1)
  svd1 <- svd(marginal.covariances$A1)#
  svd2 <- svd(marginal.covariances$A2)#
  B <- diagonal_band(Data,marginal.covariances$A1,marginal.covariances$A2) #
  lambda = svd1$d
  gamma = svd2$d
  shift.stat <- matrix(NA, l1, l2)
  for (rowi in 1:l1) {
    for (coli in 1:l2) {
      tmp <- numeric(1)
      for (i in 1:N) {
        tmp <- tmp + (as.double(t(svd1$u[, rowi]) %*%
                                  (Data[i, , ]) %*% svd2$u[, coli]))^2
      }
      tmp <- tmp/N
      shift.stat[rowi, coli] <- tmp - lambda[rowi] * gamma[coli]
      pom <- svd1$u[, rowi] %*% t(svd2$u[, coli]) #
      bandred <- sum( B*pom^2 ) #
      shift.stat[rowi, coli] <- shift.stat[rowi, coli] - bandred #
    }
  }
  shift.stat = sqrt(N) * shift.stat
  sigma.left = NULL
  sigma.right = NULL
  if (with.asymptotic.variances) {
    sigma.left <- sqrt(2) * outer(lambda[1:l1], lambda[1:l1]) *
      (diag(l1) * (sum(lambda)^2) + matrix(sum(lambda^2),
                                           l1, l1) - sum(lambda) * outer(lambda[1:l1], lambda[1:l1],
                                                                         "+"))/(sum(lambda) * sum(gamma))
    sigma.right <- sqrt(2) * outer(gamma[1:l2], gamma[1:l2]) *
      (diag(l2) * (sum(gamma)^2) + matrix(sum(gamma^2),
                                          l2, l2) - sum(gamma) * outer(gamma[1:l2], gamma[1:l2],
                                                                       "+"))/(sum(lambda) * sum(gamma))
  }
  ans = list(T.N = shift.stat, sigma.left = sigma.left, sigma.right = sigma.right)
  return(ans)
}
