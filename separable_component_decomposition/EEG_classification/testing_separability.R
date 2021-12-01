# install.packages("covsep")
library(covsep)

load("aData1.RData")
empirical_bootstrap_test(aData1,3,3)

load("aData2.RData")
empirical_bootstrap_test(aData2,3,3)

load("aData3.RData")
empirical_bootstrap_test(aData3,3,3)