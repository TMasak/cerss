setwd("C:/Users/Tomas/Documents/Skola/EPFL/Tex/PT_paper/JASA_resubmission/codes")
source("library.R")
library(ggplot2)
library(gridExtra)
library(reshape2)
library(scatterplot3d)
library(locfit)
load("mortality_data.RData")

plot_raw <- function(n){
  K1 <- dim(data)[2]; K2 <- dim(data)[3]
  s3d <- scatterplot3d(1:10, 1:10, (1:10)/10,box=F,type="n",grid=F,
                       x.ticklabs=c(1963,1973,1983,1993,2003,2013),y.ticklabs=c("60","70","80","90",""),
                       xlim=c(1,50),ylim=c(1,40),zlim = c(0,max(data[n,,])),
                       xlab="year",ylab="",zlab="mortality rate")
  text(6.7,0.2,"age",srt=22)
  text(7.45,1.8,"100",cex=0.8)
  s3d$points3d(c(0,0),c(0,K2),c(0.1,0.1),type="l",lty=3,col="gray60")
  s3d$points3d(c(0,0),c(0,K2),c(0.3,0.3),type="l",lty=3,col="gray60")
  s3d$points3d(c(0,0),c(0,K2),c(0.5,0.5),type="l",lty=3,col="gray60")
  s3d$points3d(c(0,0),c(0,K2),c(0.7,0.7),type="l",lty=3,col="gray60")
  s3d$points3d(c(0,K1),c(K2,K2),c(0.1,0.1),type="l",lty=3,col="gray60")
  s3d$points3d(c(0,K1),c(K2,K2),c(0.3,0.3),type="l",lty=3,col="gray60")
  s3d$points3d(c(0,K1),c(K2,K2),c(0.5,0.5),type="l",lty=3,col="gray60")
  s3d$points3d(c(0,K1),c(K2,K2),c(0.7,0.7),type="l",lty=3,col="gray60")
  s3d$points3d(c(0,0),c(K2,K2),c(0,0.7),type="l",lty=3,col="gray60")
  s3d$points3d(c(K1,K1),c(K2,K2),c(0,0.7),type="l",lty=3,col="gray60")
  for(i in 1:K1) s3d$points3d(rep(i, K2),1:K2, data[n,i,], type= "l")
  for(i in 1:K2) s3d$points3d(1:K1,rep(i, K1), data[n,,i], type= "l")
}

plot_raw(7)  # 7 for the Czech Republic
plot_raw(29) # 29 for Switzerland

### Cross-Validation

CV <- CV_nonstationary(data,10,5,0) # runs for some time, because when searching for d>1 without stationarity, we 
                                    # are outside of the computational limits -- complexity O(K^4)
op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(0,5),ylim=c(0.65,0.7),xaxt="n",yaxt="n",xlab="d",
     ylab="CV objective value")
axis(1, 0:5,0:5,lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
points(0:5,CV[1,]/CV[2,],type="b") # save as 5x4 inch

### Estimation

Est <- estimate_separable(data,1)
B <- diagonal_band(data,Est$A1,Est$A2)

plotB <- ggplot(data = melt(log(B+0.01)), mapping = aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + theme_minimal() + 
  theme(text=element_text(size=14, family="serif"), plot.title = element_text(hjust = 0.5, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
  scale_x_continuous(breaks=c(1,10,20,30,40,50),labels=c(1964,1973,1983,1993,2003,2013),expand = c(0,0)) + 
  scale_y_continuous(breaks=c(1,10,20,30,40),labels=c(60,70,80,90,100),expand = c(0,0)) +
  xlab(label = "year") + ylab(label = "age") +
  scale_fill_gradient(low='mediumblue', high ='yellow',
                      guide = guide_colorbar(barwidth=1, barheight = 17)) + labs(fill="")
plotB # save as 4.65x3.7 inch

### Eigendecompositions

EIG1_shift <- eigen(Est$A1)
EIG1_shift$values <- pmax(EIG1_shift$values,0)
EIG2_shift <- eigen(Est$A2)
EIG2_shift$values <- pmax(EIG2_shift$values,0)
Est <- estimate_separable(data,0)
EIG1 <- eigen(Est$A1)
EIG2 <- eigen(Est$A2)

sum(EIG1$values[1:16])/sum(EIG1$values) # 16 eigenfunctions needed to capture 90% of non-shifted A1
sum(EIG2$values[1:4])/sum(EIG2$values)  # 4 eigenfunctions needed to capture 90% of non-shifted A1
sum(EIG1_shift$values[1:4])/sum(EIG1_shift$values) # 4 eigenfunctions needed to capture 90% of shifted A1
sum(EIG2_shift$values[1:2])/sum(EIG2_shift$values) # 2 eigenfunctions needed to capture 90% of non-shifted A1

sum(EIG1_shift$values[1:2])/sum(EIG1_shift$values) # 2 eigenfunctions of shifted A1 capture 83%

op <- par(mfrow=c(1,1),mar = c(3.2, 2, 1.6, 0.4), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,50),ylim=c(-1,1),xaxt="n",yaxt="n",xlab="year",ylab="")
axis(1,0:5*10,1964+0:5*10,lwd=0.7,mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
points(-EIG1_shift$vectors[,1],type="l")
points(-EIG1$vectors[,1],type="l",lty=2)
legend("bottomleft",legend=c("SPT","PT"),lty=c(1,2))

op <- par(mfrow=c(1,1),mar = c(3.2, 2, 1.6, 0.4), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,50),ylim=c(-1,1),xaxt="n",yaxt="n",xlab="year",ylab="")
axis(1,0:5*10,1964+0:5*10,lwd=0.7,mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
points(-EIG1_shift$vectors[,2],type="l")
points(-EIG1$vectors[,2],type="l",lty=2)
legend("bottomleft",legend=c("SPT","PT"),lty=c(1,2))

op <- par(mfrow=c(1,1),mar = c(3.2, 2, 1.6, 0.4), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,40),ylim=c(-1,1),xaxt="n",yaxt="n",xlab="age",ylab="")
axis(1,0:4*10,60+0:4*10,lwd=0.7,mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
points(-EIG2_shift$vectors[,1],type="l")
points(-EIG2$vectors[,1],type="l",lty=2)
legend("bottomleft",legend=c("SPT","PT"),lty=c(1,2))

op <- par(mfrow=c(1,1),mar = c(3.2, 2, 1.6, 0.4), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,40),ylim=c(-1,1),xaxt="n",yaxt="n",xlab="age",ylab="")
axis(1,0:4*10,60+0:4*10,lwd=0.7,mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
points(-EIG2_shift$vectors[,2],type="l")
points(EIG2$vectors[,2],type="l",lty=2) # save as 4x4
legend("bottomleft",legend=c("SPT","PT"),lty=c(1,2))


