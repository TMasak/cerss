setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("library_PIP.R")

papcol <- c("forestgreen","royalblue","orange","magenta","red","saddlebrown")
pappch <- c(0,4,2,19,18,5)
paplty <- rep(1,6) #c(4,2,3,1,5,6)
legenda <- function(altitude,dash,symb,barva,nazev){
  points(c(start,start+width),c(altitude,altitude),type="l",lty=dash,lwd=2,col=barva)
  points(c(start+3/2*width),c(altitude),pch=symb,lwd=2,col=barva)
  points(c(start+2*width,start+3*width),c(altitude,altitude),type="l",lty=dash,lwd=2,col=barva)
  text(start+space,altitude,nazev,adj=0)
}

################
### Gneiting ###
################

K1 <- 50
K2 <- 50
beta <- 0.7
C <- get_gneiting_cov(K1,K2,beta)

C_mat <- tensor2matrix(C)
C_mat_perm <- VanLoansPerm(C_mat,K1,K2)
EIG <- eigen(C_mat)
SVD <- svd(C_mat_perm)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(EIG$values[1:1000],type="h",ylab="Magnitude")
plot(SVD$d[1:9],type="h",ylab="Magnitude",lwd=2,ylim=c(0,200),xaxt="n")
axis(1,at=seq(1,9,by=2),labels=seq(1,9,by=2))

load("data_gneiting_estimation_emp.RData")
Results <- Data[[1]]
for(i in 2:10){
  Akt <- Data[[i]]
  Results[,,i] <- Akt[,,i]
}
Results <- Results[,,1:10]
meanResults <- apply(Results,c(1,2),mean)
EMP <- meanResults[,2]

load("data_gneiting_estimation.RData")
Results <- Data[[1]]
for(i in 2:10){
  Akt <- Data[[i]]
  Results[,,i] <- Akt[,,i]
}
Results <- Results[,,1:10]
meanResults <- apply(Results,c(1,2),mean)
CV2picks <- Results[,16,]
CV2errors <- CV2picks
for(i in 1:9){
  for(j in 1:10){
    CV2errors[i,j] <- Results[i,CV2picks[i,j],j]
  }
}
meanCV2errors <- rowMeans(CV2errors) 

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.25), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(log(2^(7:15),2),meanResults[,1],type="b",ylim=c(0,0.3),xaxt="n",xlab="N",
     ylab="relative error",col="white")
# abline(v=log(2^c(8,11),2),col="gray60")
# lines(x=c(log(2^14,2),log(2^14,2)),y=c(-1,0.187),col="gray60")
axis(1,at=log(2^(7:15),2),labels = 2^(7:15))
abline(h=sqrt(sum(SVD$d[2:2500]^2))/sqrt(sum(SVD$d^2)),lty=2,lwd=1.5,col=papcol[1])
abline(h=sqrt(sum(SVD$d[3:2500]^2))/sqrt(sum(SVD$d^2)),lty=2,lwd=1.5,col=papcol[2])
abline(h=sqrt(sum(SVD$d[4:2500]^2))/sqrt(sum(SVD$d^2)),lty=2,lwd=1.5,col=papcol[3])
# abline(h=sqrt(sum(SVD$d[5:2500]^2))/sqrt(sum(SVD$d^2)),col="forestgreen",lty=5)
abline(h=0,lty=2,lwd=1.5,col=papcol[5])
points(log(2^(7:15),2),meanResults[,1],type="b",pch=pappch[1],lwd=2,lty=paplty[1],col=papcol[1])
points(log(2^(7:15),2),meanResults[,2],type="b",pch=pappch[2],lwd=2,lty=paplty[2],col=papcol[2])
points(log(2^(7:15),2),meanResults[,3],type="b",pch=pappch[3],lwd=2,lty=paplty[3],col=papcol[3])
points(log(2^(7:15),2),meanCV2errors,type="b",ylim=c(0,0.3),pch=pappch[4],lwd=2,lty=paplty[4],col=papcol[4])
points(log(2^(7:15),2),EMP,type="b",pch=pappch[5],lwd=2,lty=paplty[5],col=papcol[5])

### manual legend
space <- 1
width <- 0.25
height <- 2
start <- 12.6
legenda(0.3,  paplty[1],pappch[1],papcol[1],"separable")
legenda(0.275,paplty[2],pappch[2],papcol[2],"2-separable")
legenda(0.25, paplty[3],pappch[3],papcol[3],"3-separable")
legenda(0.225,paplty[4],pappch[4],papcol[4],"cross-validated")
legenda(0.2,  paplty[5],pappch[5],papcol[5],"empirical")
lines(12.5+c(0,2.8,2.8,0,0),0.187+c(0,0,0.124,0.124,0)) # rectangle
# save as 7.7x4.83

### prediction
load("data_gneiting_prediction_new.RData")
Results <- array(0,c(dim(Data[[1]])[1],dim(Data[[1]])[2],25))
for(i in 1:25){
  Results[,,i] <- Data[[i]]
}
meanResults <- apply(Results,c(1,2),mean)
round(t(meanResults[,c(16,17,18,20)]),digits=3)

######################
### rancov fourier ###
###################### 

load("data_fourier_Chen_99.Rdata")
Data1 <- Data
load("data_fourier.RData")
for(i in 1:25){
  for(j in 1:dim(Data[[1]])[1]){
    Data[[i]][j,4] <- Data[[i]][j,Data[[i]][j,13]]
    Data[[i]][j,19] <- Data[[i]][j,Data[[i]][j,13]+15]
    Data[[i]][j,24] <- Data[[i]][j,Data[[i]][j,13]+20]
    Data[[i]][j,14] <- Data1[[i]][j,14]
    Data[[i]][j,15] <- Data1[[i]][j,25]
  }
}
Results <- Data[[1]]
for(i in 2:25){
  Results <- Results + Data[[i]]
}
Results <- Results/25
x <- 2^(3:1)
Bias <- rep(0,3)
for(i in 1:3) Bias[i] <- sqrt(sum(x[(i+1):3]^2)/sum(x^2))

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.25), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(log(2^(7:15),2),Results[,1],type="b",ylim=c(0,1),xaxt="n",xlab="N",
     ylab="relative error",col="white")
axis(1,at=log(2^(7:15),2),labels = 2^(4:12))
points(log(2^(7:15),2),Results[,1],type="b",pch=pappch[1],lwd=2,lty=paplty[1],col=papcol[1])
points(log(2^(7:15),2),Results[,2],type="b",pch=pappch[2],lwd=2,lty=paplty[2],col=papcol[2])
points(log(2^(7:15),2),Results[,3],type="b",pch=pappch[3],lwd=2,lty=paplty[3],col=papcol[3])
points(log(2^(7:15),2),Results[,4],type="b",pch=pappch[4],lwd=2,lty=paplty[4],col=papcol[4])
points(log(2^(7:15),2),Results[,5],type="b",pch=pappch[5],lwd=2,lty=paplty[5],col=papcol[5])
points(log(2^(7:15),2),Results[,14],type="b",pch=pappch[6],lwd=2,lty=paplty[6],col=papcol[6]) # Chen
abline(h=Bias[1],lty=2,lwd=1.5,col=papcol[1])
abline(h=Bias[2],lty=2,lwd=1.5,col=papcol[2])
abline(h=0,lty=2,lwd=1.5,col=papcol[3])
# abline(h=0,lty=2,lwd=1.5,col=papcol[5])
space <- 1
width <- 0.25
height <- 2
start <- 12.4
legenda(1,  paplty[1],pappch[1],papcol[1],"separable")
legenda(0.92,paplty[2],pappch[2],papcol[2],"2-separable")
legenda(0.84, paplty[3],pappch[3],papcol[3],"3-separable")
legenda(0.76,paplty[4],pappch[4],papcol[4],"cross-validated")
legenda(0.68,  paplty[5],pappch[5],papcol[5],"empirical")
legenda(0.6,  paplty[6],pappch[6],papcol[6],"weakly separable")
lines(12.3+c(0,2.95,2.95,0,0),0.56+c(0,0,0.48,0.48,0)) # rectangle

plot(log(2^(7:15),2),Results[,1],type="b",ylim=c(0,1),xaxt="n",xlab="N",
     ylab="relative error",col="white")
axis(1,at=log(2^(7:15),2),labels = 2^(4:12))
abline(h=c(0.2,0.4,0.6,0.8,1),col="gray60")
points(log(2^(7:15),2),Results[,16+5],type="b",pch=pappch[1],lwd=2,lty=paplty[1],col=papcol[1])
points(log(2^(7:15),2),Results[,17+5],type="b",pch=pappch[2],lwd=2,lty=paplty[2],col=papcol[2])
points(log(2^(7:15),2),Results[,18+5],type="b",pch=pappch[3],lwd=2,lty=paplty[3],col=papcol[3])
points(log(2^(7:15),2),Results[,19+5],type="b",pch=pappch[4],lwd=2,lty=paplty[4],col=papcol[4])
points(log(2^(7:15),2),Results[,20+5],type="b",pch=pappch[5],lwd=2,lty=paplty[5],col=papcol[5])
points(log(2^(7:15),2),Results[,15],type="b",pch=pappch[6],lwd=2,lty=paplty[6],col=papcol[6]) # Chen
space <- 1
width <- 0.25
height <- 2
start <- 6.8
legenda(0.16,  paplty[1],pappch[1],papcol[1],"separable")
legenda(0.08,paplty[2],pappch[2],papcol[2],"2-separable")
legenda(0   , paplty[3],pappch[3],papcol[3],"3-separable")
start <- 9.3
legenda(0.16,paplty[4],pappch[4],papcol[4],"cross-validated")
legenda(0.08,  paplty[5],pappch[5],papcol[5],"empirical")
legenda(0   , paplty[6],pappch[6],papcol[6],"weakly separable")
lines(6.68+c(0,5.55,5.55,0,0),-0.04+c(0,0,0.24,0.24,0)) # rectangle

### runtimes

load("data_fourier_times.RData")
op <- par(mfrow=c(1,1),mar = c(3.2, 4, 1.6, 0.25), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(seq(50,350,by=30),seq(50,350,by=30),type="b",ylim=c(0,7),xaxt="n",yaxt="n",xlab="K",
     ylab="",col="white")
axis(1,at=seq(50,350,by=30),labels = seq(50,350,by=30))
axis(2,at=c(0,0.6931472,1.791759,2.772589,3.931826,5.303305,7.09091),labels=c(0,1,5,15,50,200,1200),las=1)
abline(h=c(0,0.6931472,1.791759,2.772589,3.931826),col="gray60")
points(c(0,255),c(5.303305,5.303305),type="l",col="gray60")
points(c(0,255),c(7.09091,7.09091),type="l",col="gray60")
abline(v=140,col="gray60")
title(ylab="runtime [sec]",line=3)
points(seq(50,350,by=30),log(1+ERR[,6]),type="b",pch=pappch[1],lwd=2,lty=paplty[1],col=papcol[1])
points(seq(50,350,by=30),log(1+ERR[,7]),type="b",pch=pappch[2],lwd=2,lty=paplty[2],col=papcol[2])
points(seq(50,350,by=30),log(1+ERR[,8]),type="b",pch=pappch[3],lwd=2,lty=paplty[3],col=papcol[3])
points(seq(50,140,by=30),log(1+25*ERR[1:4,10]),type="b",pch=pappch[5],lwd=2,lty=paplty[5],col=papcol[5])
space <- 1*50
width <- 0.25*50
height <- 2*50
start <- 260
legenda(7,  paplty[1],pappch[1],papcol[1],"separable")
legenda(6.5,paplty[2],pappch[2],papcol[2],"2-separable")
legenda(6, paplty[3],pappch[3],papcol[3],"3-separable")
legenda(5.5,  paplty[5],pappch[5],papcol[5],"empirical")
lines(255+c(0,106.5,106.5,0,0),5.28+c(0,0,1.97,1.97,0)) # rectangle

#########################
### Random Covariance ###
#########################

load("data_rancov.RData")
Results <- Data[[1]]
for(i in 2:10){
  Results <- Results + Data[[i]]
}
Results <- Results/10

Bias <- rep(0,9)
for(i in 1:9){
  x <- (i+1)^(4:1)
  Bias[i] <- sqrt(x[3]^2+x[4]^2)/sqrt(sum(x^2))
}

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(log(2^(7:13),2),Results[,1,2]-Bias[1],type="b",ylim=c(0,0.3),xaxt="n",xlab="N",
     ylab="bias-free relative error")
axis(1,at=log(2^(7:13),2),labels = 2^(7:13))
points(log(2^(7:13),2),Results[,2,2]-Bias[2],type="b",col="blue",pch=0)
points(log(2^(7:13),2),Results[,3,2]-Bias[3],type="b",col="darkviolet",pch=2)
points(log(2^(7:13),2),Results[,4,2]-Bias[4],type="b",col="red",pch=3)
# points(log(2^(7:13),2),Results[,5,2]-Bias[5],type="b",col="orange")
points(log(2^(7:13),2),Results[,6,2]-Bias[6],type="b",col="forestgreen",pch=4)
# points(Results[,7,2]-Bias[7],type="b",col="brown")
legend("topright",legend=c(expression(paste(alpha,"=2")),
                           expression(paste(alpha,"=3")),
                           expression(paste(alpha,"=4")),
                           expression(paste(alpha,"=5")),
                           expression(paste(alpha,"=6      "))),
       pch=c(1,0,2,3,4),lwd=1.5,col=c("black","blue","darkviolet","red","forestgreen"), y.intersp=1.5)

Bias <- rep(0,9)
for(i in 1:9){
  x <- (i+1)^(4:1)
  Bias[i] <- x[4]/sqrt(sum(x^2))
}

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(log(2^(7:13),2),Results[,1,3]-Bias[1],type="b",ylim=c(0,0.3),xaxt="n",xlab="N",
     ylab="bias-free relative error")
axis(1,at=log(2^(7:13),2),labels = 2^(7:13))
points(log(2^(7:13),2),Results[,2,3]-Bias[2],type="b",col="blue",0)
points(log(2^(7:13),2),Results[,3,3]-Bias[3],type="b",col="darkviolet",2)
points(log(2^(7:13),2),Results[,4,3]-Bias[4],type="b",col="red",3)
# points(Results[,5,3]-Bias[5],type="b",col="orange")
points(log(2^(7:13),2),Results[,6,3]-Bias[6],type="b",col="forestgreen",4)
# points(Results[,7,3]-Bias[7],type="b",col="brown")
legend("topright",legend=c(expression(paste(alpha,"=2")),
                           expression(paste(alpha,"=3")),
                           expression(paste(alpha,"=4")),
                           expression(paste(alpha,"=5")),
                           expression(paste(alpha,"=6      "))),
       pch=c(1,0,2,3,4),lwd=1.5,col=c("black","blue","darkviolet","red","forestgreen"), y.intersp=1.5)

#######################
### Inverse Problem ###
#######################

load("data_iteration2.RData")

Results <- array(0,c(100,7,3,2))
for(i in 1:100){
  Results[i,,,] <- Data[[i]]
}
Results[,,2,] <- Results[,,1,]

load("data_iteration1.RData")
for(i in 1:100){
  Akt <- Data[[i]]
  Results[i,,1,] <- Akt[,1,]
}

max(Results[,,,1]) # check that inverse algorithm never failed
Results <- Results[,,,2]
meanResults <- apply(Results,c(2,3),mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(1:7,meanResults[,1],type="b",ylim=c(0,50),xaxt="n",xlab="grid size", ylab="no. of iterations",
     lwd=1.5)
axis(1,at=1:7,labels = (1:7)*30)
points(meanResults[,2],type="b",col="blue",pch=0,lwd=1.5)
points(meanResults[,3],type="b",col="red",pch=2,lwd=1.5)
legend("topright",legend=c(expression(paste(kappa,"=10")),
                           expression(paste(kappa,"=100")),
                           expression(paste(kappa,"=1000      "))),
       pch=c(1,0,2),lwd=1.5,col=c("black","blue","red"), y.intersp=1.5)

load("data_condition.RData")
Results <- apply(Data,2:5,mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 0.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(Results[1,1,,2],Results[1,1,,1],type="b",ylim=c(0,1),xlim=c(0,60),
     xlab="no. of iterations",ylab=expression(sigma[1]),
     col="brown4",yaxt="n") # rank=7, eps=0.5
axis(2,at=seq(0.15,0.95,by=0.2),labels=seq(0.15,0.95,by=0.2))
points(Results[2,1,1:8,2],Results[2,1,1:8,1],type="b",pch=4,col="brown4") # rank=5
points(Results[3,1,1:7,2],Results[3,1,1:7,1],type="b",pch=2,col="brown4") # rank=3
points(Results[1,2,,2],Results[1,2,,1],type="b",col="blue") # rank=7, eps=0.2
points(Results[2,2,1:8,2],Results[2,2,1:8,1],type="b",col="blue",pch=4) # rank=5
points(Results[3,2,1:7,2],Results[3,2,1:7,1],type="b",col="blue",pch=2) # rank=3
points(Results[1,3,,2],Results[1,3,,1],type="b",col="red") # rank=7, eps=0.1
points(Results[2,3,1:8,2],Results[2,3,1:8,1],type="b",col="red",pch=4) # rank=5
points(Results[3,3,1:7,2],Results[3,3,1:7,1],type="b",col="red",pch=2) # rank=3
points(Results[1,4,,2],Results[1,4,,1],type="b",col="forestgreen") # rank=7, eps=0.05
points(Results[2,4,1:8,2],Results[2,4,1:8,1],type="b",col="forestgreen",pch=4) # rank=5
points(Results[3,4,1:7,2],Results[3,4,1:7,1],type="b",col="forestgreen",pch=2) # rank=3
points(Results[1,5,,2],Results[1,5,,1],type="b",col="magenta") # rank=7, eps=0.02
points(Results[2,5,1:8,2],Results[2,5,1:8,1],type="b",col="magenta",pch=4) # rank=5
points(Results[3,5,1:7,2],Results[3,5,1:7,1],type="b",col="magenta",pch=2) # rank=3
points(Results[1,6,,2],Results[1,6,,1],type="b",col="darkorange") # rank=7, eps=0.01
points(Results[2,6,1:8,2],Results[2,6,1:8,1],type="b",col="darkorange",pch=4) # rank=5
points(Results[3,6,1:7,2],Results[3,6,1:7,1],type="b",col="darkorange",pch=2) # rank=3
legend("left",legend=c("R=3","R=5","R=7   "),pch=c(2,4,1), y.intersp=1.5)
legend("bottomleft",legend=c(expression(paste(epsilon,"=0.5")),
                             expression(paste(epsilon,"=0.2")),
                             expression(paste(epsilon,"=0.1")),
                             expression(paste(epsilon,"=0.05")),
                             expression(paste(epsilon,"=0.02")),
                             expression(paste(epsilon,"=0.01      "))),
       col=c("brown4","blue","red","forestgreen","magenta","darkorange"),lty=1, y.intersp=1.5)


### Weakly Separable Model - comparison of different cutoffs for different eigenspace dimensions (for review only)

load("data_fourier_Chen.Rdata")
Data1 <- Data
load("data_fourier_Chen_3.Rdata")
Data2 <- Data
load("data_fourier_Chen_50.Rdata")
Data3 <- Data
load("data_fourier_Chen_90.Rdata")
Data4 <- Data
load("data_fourier_Chen_99.Rdata")
for(i in 1:25){
  for(j in 1:dim(Data[[1]])[1]){
    Data[[i]][j,10] <- Data1[[i]][j,14]
    Data[[i]][j,11] <- Data2[[i]][j,14]
    Data[[i]][j,12] <- Data3[[i]][j,14]
    Data[[i]][j,13] <- Data4[[i]][j,14]
  }
}
Results <- Data[[1]]
for(i in 2:25){
  Results <- Results + Data[[i]]
}
Results <- Results/25

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 0.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(log(2^(7:15),2),Results[,10],ylim=c(0,1),type="b",lwd=2,pch=pappch[6],lty=paplty[6],col=papcol[6],xaxt="no",
     xlab="N")
axis(1,at=log(2^(7:15),2),labels = 2^(4:12))
abline(h=seq(0,1,by=0.2),col="gray60")
points(log(2^(7:15),2),Results[,11],ylim=c(0,1),type="b",lwd=1.5,pch=4,col="red")
# points(log(2^(7:15),2),Results[,12],ylim=c(0,1),type="b",lwd=1.5,pch=0)
points(log(2^(7:15),2),Results[,13],ylim=c(0,1),type="b",lwd=1.5,pch=3,col="blue")
points(log(2^(7:15),2),Results[,14],ylim=c(0,1),type="b",lwd=1.5,pch=6,col="forestgreen")
points(log(2^(7:15),2),Results[,10],ylim=c(0,1),type="b",lwd=2,pch=pappch[6],lty=paplty[6],col=papcol[6])
legend("right",legend=c("I,J=3","90 %","95 % (used)", "99 %"),pch=c(4,3,pappch[6],6),
       col=c("red","blue",papcol[6],"forestgreen"),pt.lwd = c(1.5,1.5,2,1.5))



