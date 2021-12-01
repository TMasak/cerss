setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("./data/surfaces.RData")
load("./data/trial_info.RData")
library(lattice)
# Data <- Data/100

alcoholic <- as.logical(c(rep(1,65),rep(0,44),rep(1,12),0))

Data <- Data[-47,,,]
Data2 <- Data2[-47,]
alcoholic <- alcoholic[-47]

cond1 <- I(Data2 == 1)
cond2 <- I(Data2 == 2)
cond3 <- I(Data2 == 3)

set.seed(517)
aData1 <- aData2 <- aData3 <- array(0,c(121,64,256))
for(i in 1:121){
  ind <- sample(which(cond1[i,]),10) # take random 10 trials for condition 1
  aData1[i,,] <- apply(Data[i,ind,,],c(2,3),mean)
  ind <- sample(which(cond2[i,]),10)
  aData2[i,,] <- apply(Data[i,ind,,],c(2,3),mean)
  ind <- sample(which(cond3[i,]),10)
  aData3[i,,] <- apply(Data[i,ind,,],c(2,3),mean)
}

save(aData1,file="aData1.RData")
save(aData2,file="aData1.RData")
save(aData3,file="aData1.RData")