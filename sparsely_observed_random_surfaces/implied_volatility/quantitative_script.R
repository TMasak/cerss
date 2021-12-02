setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load("grid20.RData")

Ratios <- array(0,c(392,5))
for(j in 1:5){
  Ratios[,j] <- c(rmse_our_abcde[[j]] / rmse_smoother_abcde[[j]],
                  rmse_4D_abcde[[j]] / rmse_smoother_abcde[[j]])
}
ratio <- c(Ratios)

load("TRout.RData")
TRout$ratio <- ratio
save(TRout,file="TRout.RData")
