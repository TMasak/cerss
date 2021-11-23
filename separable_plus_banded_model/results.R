setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('library.R')
load("Norms.RData")
Norms[1,1:4] <- 1
setwd("./data")
library(ggplot2)
library(reshape2)

###################################
### Brownian and Legendre plots ###
###################################

load("brownian_ar_buf.RData")
Dat <- array(0,c(25,11,7))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(0,19),ylim=c(0,1),xaxt="n",yaxt="n",xlab="d",
     ylab="relative estimation error")
axis(1, c(-0.1,1.1,seq(3,19,2)),c(0,seq(1,19,2)),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
abline(v=9,lty=3)
for(j in c(1,2,3,5,7)) points(c(0,seq(1,19,2)),meanDat[,j],type="b",
                     col=papcol(j),pch=pappch(j),lwd=1.5)
points(c(0,seq(1,19,2)),1-Norms[,1],type="b", col=papcol(6),pch=pappch(6),lwd=1.5)
# save as 5x4 inch

load("brownian_tri_buf.RData")
Dat <- array(0,c(25,11,7))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(0,19),ylim=c(0,1),xaxt="n",yaxt="n",xlab="d",
     ylab="relative estimation error")
axis(1, c(-0.1,1.1,seq(3,19,2)),c(0,seq(1,19,2)),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
abline(v=9,lty=3)
for(j in c(1,2,3,5,7)) points(c(0,seq(1,19,2)),meanDat[,j],type="b",
                                col=papcol(j),pch=pappch(j),lwd=1.5)
points(c(0,seq(1,19,2)),1-Norms[,2],type="b", col=papcol(6),pch=pappch(6),lwd=1.5)
# save as 5x4 inch

load("legendre_ar_buf.RData")
Dat <- array(0,c(25,11,7))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)
# rel <- array(0,c(10,6))
# for(j in 1:5) rel[,j] <- meanDat[,6]/meanDat[,j]

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(0,19),ylim=c(0,1),xaxt="n",yaxt="n",xlab="d",
     ylab="relative estimation error")
axis(1, c(-0.1,1.1,seq(3,19,2)),c(0,seq(1,19,2)),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
abline(v=9,lty=3)
for(j in c(1,2,3,5,7)) points(c(0,seq(1,19,2)),meanDat[,j],type="b",
                                col=papcol(j),pch=pappch(j),lwd=1.5)
points(c(0,seq(1,19,2)),1-Norms[,3],type="b", col=papcol(6),pch=pappch(6),lwd=1.5)
# save as 5x4 inch

load("legendre_tri_buf.RData")
Dat <- array(0,c(25,11,7))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)
# rel <- array(0,c(10,6))
# for(j in 1:5) rel[,j] <- meanDat[,6]/meanDat[,j]

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(0,19),ylim=c(0,1),xaxt="n",yaxt="n",xlab="d",
     ylab="relative estimation error")
axis(1, c(-0.1,1.1,seq(3,19,2)),c(0,seq(1,19,2)),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
abline(v=9,lty=3)
for(j in c(1,2,3,5,7)) points(c(0,seq(1,19,2)),meanDat[,j],type="b",
                                col=papcol(j),pch=pappch(j),lwd=1.5)
points(c(0,seq(1,19,2)),1-Norms[,4],type="b", col=papcol(6),pch=pappch(6),lwd=1.5)
# save as 5x4 inch

#############################################################################

load("brownian_ar_tau.RData")
Dat <- array(0,c(25,9,7))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,9),ylim=c(0,1),xaxt="n",yaxt="n",xlab=expression(tau),
     ylab="relative estimation error")
axis(1, 1:9,c(1,2,4,8,16,32,64,128,256),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
abline(v=log(3,2)+1,lty=3)
for(j in c(1,2,3,5,7)) points(1:9,meanDat[,j],type="b",col=papcol(j),pch=pappch(j),lwd=1.5)
points(1:9,1-Norms[1:9,5],type="b", col=papcol(6),pch=pappch(6),lwd=1.5)
legend("topright",legend=c("ECE","NKP","PT","SPT","SPT-CV","bias"),col=papcol(c(7,5,2,1,3,6)),pch=pappch(c(7,5,2,1,3,6)),pt.lwd=1.5)

load("brownian_tri_tau.RData")
Dat <- array(0,c(25,9,7))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,9),ylim=c(0,1),xaxt="n",yaxt="n",xlab=expression(tau),
     ylab="relative estimation error")
axis(1, 1:9,c(1,2,4,8,16,32,64,128,256),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
abline(v=log(3,2)+1,lty=3)
for(j in c(1,2,3,5,7)) points(1:9,meanDat[,j],type="b",col=papcol(j),pch=pappch(j),lwd=1.5)
points(1:9,1-Norms[1:9,6],type="b", col=papcol(6),pch=pappch(6),lwd=1.5)
legend("topright",legend=c("ECE","NKP","PT","SPT","SPT-CV","bias"),col=papcol(c(7,5,2,1,3,6)),pch=pappch(c(7,5,2,1,3,6)),pt.lwd=1.5)

load("legendre_ar_tau.RData")
Dat <- array(0,c(25,9,7))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,9),ylim=c(0,1),xaxt="n",yaxt="n",xlab=expression(tau),
     ylab="relative estimation error")
axis(1, 1:9,c(1,2,4,8,16,32,64,128,256),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
abline(v=log(3,2)+1,lty=3)
for(j in c(1,2,3,5,7)) points(1:9,meanDat[,j],type="b",col=papcol(j),pch=pappch(j),lwd=1.5)
points(1:9,1-Norms[1:9,7],type="b", col=papcol(6),pch=pappch(6),lwd=1.5)
legend("topright",legend=c("ECE","NKP","PT","SPT","SPT-CV","bias"),col=papcol(c(7,5,2,1,3,6)),pch=pappch(c(7,5,2,1,3,6)),pt.lwd=1.5)

load("legendre_tri_tau.RData")
Dat <- array(0,c(25,9,7))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,9),ylim=c(0,1),xaxt="n",yaxt="n",xlab=expression(tau),
     ylab="relative estimation error")
axis(1, 1:9,c(1,2,4,8,16,32,64,128,256),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
abline(v=log(3,2)+1,lty=3)
for(j in c(1,2,3,5,7)) points(1:9,meanDat[,j],type="b",col=papcol(j),pch=pappch(j),lwd=1.5)
points(1:9,1-Norms[1:9,8],type="b", col=papcol(6),pch=pappch(6),lwd=1.5)
legend("topright",legend=c("ECE","NKP","PT","SPT","SPT-CV","bias"),col=papcol(c(7,5,2,1,3,6)),pch=pappch(c(7,5,2,1,3,6)),pt.lwd=1.5)

#############################################################################

load("brownian_ar_N.RData")
Dat <- array(0,c(25,9,7))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,7),ylim=c(0,1),xaxt="n",yaxt="n",xlab="N",
     ylab="relative estimation error")
axis(1, 1:7,2^(5:11),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
abline(v=log(300,2)-4,lty=3)
for(j in c(1,2,3,5,7)) points(1:7,meanDat[1:7,j],type="b",col=papcol(j),pch=pappch(j),lwd=1.5)
points(1:7,1-Norms[1:7,9],type="b", col=papcol(6),pch=pappch(6),lwd=1.5)
# save as 5x4 inch

load("brownian_tri_N.RData")
Dat <- array(0,c(25,9,7))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,7),ylim=c(0,1),xaxt="n",yaxt="n",xlab=expression(tau),
     ylab="relative estimation error")
axis(1, 1:7,2^(5:11),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
abline(v=log(300,2)-4,lty=3)
for(j in c(1,2,3,5,7)) points(1:7,meanDat[1:7,j],type="b",col=papcol(j),pch=pappch(j),lwd=1.5)
points(1:7,1-Norms[1:7,10],type="b", col=papcol(6),pch=pappch(6),lwd=1.5)
# save as 5x4 inch

load("legendre_ar_N.RData")
Dat <- array(0,c(25,9,7))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,7),ylim=c(0,1),xaxt="n",yaxt="n",xlab="N",
     ylab="relative estimation error")
axis(1, 1:7,2^(5:11),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
abline(v=log(300,2)-4,lty=3)
for(j in c(1,2,3,5,7)) points(1:7,meanDat[1:7,j],type="b",col=papcol(j),pch=pappch(j),lwd=1.5)
points(1:7,1-Norms[1:7,11],type="b", col=papcol(6),pch=pappch(6),lwd=1.5)
# save as 5x4 inch

load("legendre_tri_N.RData")
Dat <- array(0,c(25,9,7))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,7),ylim=c(0,1),xaxt="n",yaxt="n",xlab="N",
     ylab="relative estimation error")
axis(1, 1:7,2^(5:11),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
abline(v=log(300,2)-4,lty=3)
for(j in c(1,2,3,5,7)) points(1:7,meanDat[1:7,j],type="b",col=papcol(j),pch=pappch(j),lwd=1.5)
points(1:7,1-Norms[1:7,12],type="b", col=papcol(6),pch=pappch(6),lwd=1.5)
# save as 5x4 inch

################
### ADI runs ###
################

load("data_asi.RData")

Dat <- array(0,c(25,10,9))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)
meanDat[,7] <- meanDat[,5];meanDat[,5] <- meanDat[,4];


op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(30,210),ylim=c(0,1),xaxt="n",yaxt="n",xlab="K",
     ylab="relative estimation error")
axis(1, seq(30,210,20),seq(30,210,20),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
# for(j in 1:6) points(seq(30,210,20),meanDat[,j],type="b",col=papcol(j),pch=pappch(j),lty=paplty(j),lwd=1.5)
for(j in c(1,2,3,5,7)) points(seq(30,210,20),meanDat[,j],type="b",col=papcol(j),pch=pappch(j),lwd=1.5)
# save as 5x4 inch

legend("topright",legend=c("ECE","PT","NKP","SPT","SPT-CV"),col=papcol(c(7,2,5,1,3)),pch=pappch(c(7,2,5,1,3)),pt.lwd=1.5)

### druhy plot

load("data_adi.RData")
Dat <- array(0,c(25,10,9))
for(j in 1:25){
  Dat[j,,] <- data[[j]]
}
meanDat <- apply(Dat,c(2,3),mean)

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(30,210),ylim=c(0,50),xaxt="n",yaxt="n",xlab="K",
     ylab="iterations")
axis(1, seq(30,210,20),seq(30,210,20),lwd=0.7,
     mar = c(2, 3, 0.5, 0) + 0.5)
axis(2,lwd=0.7)
grid(nx=NA,ny=NULL,lwd=0.7,col="gray80",lty=1)
points(seq(30,210,20),meanDat[,8],type="b",pch=2,lwd=1.5)
points(seq(30,210,20),meanDat[,9],type="b",pch=3,lwd=1.5)
legend("right",legend=c("ADI","PCG"),pch=c(2,3))

########################
### Simulation Setup ###
########################

# shifted_traces <- function(C){
#   K <- dim(C)[1]
#   RES <- rep(0,K)
#   for(delta in 0:(K-1)){
#     indx <- 1:(K-delta) 
#     W1 <- matrix(0,K,K)
#     W1[cbind(indx+delta,indx)] <- W1[cbind(indx,indx+delta)] <- 1
#     RES[delta+1] <- sum(W1*C)/2
#   }
#   RES[1] <- 2*RES[1]
#   return(RES)
# }

Dat <- create_data(seed=17,buf=4,A1="brownian",A2="legendre",mask="triangular",tau=0.5)
discrete <- 30

X <- Dat$A1
# X <- cov2cor(X)
max <- max(X)
X <- round(X/max*discrete)*max/discrete
X <- melt(X)
mine.heatmap <- ggplot(data = X, mapping = aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + theme_minimal() +
  theme(text=element_text(size=14, family="serif"), plot.title = element_text(hjust = 0.5, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
  scale_x_continuous(breaks=c(1,20,40,60,80,100),labels=c(0,0.2,0.4,0.6,0.8,1),expand = c(0,0)) + 
  scale_y_continuous(breaks=c(1,20,40,60,80,100),labels=c(0,0.2,0.4,0.6,0.8,1),expand = c(0,0)) +
  coord_cartesian(xlim = c(1,100), ylim=c(1,100)) +
  xlab(label = "t (or s)") + ylab(label = "t' (or s')") +
  scale_fill_gradient(low='mediumblue', high ='yellow',
                      guide = guide_colorbar(barwidth=1, barheight = 18.5))
mine.heatmap # save as 5x4 inch

# plot(shifted_traces(X)^2,type="b")

X <- Dat$A2
# X <- cov2cor(X)
max <- max(X)
X <- round(X/max*discrete)*max/discrete
X <- melt(X)
mine.heatmap <- ggplot(data = X, mapping = aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + theme_minimal() +
  theme(text=element_text(size=14, family="serif"), plot.title = element_text(hjust = 0.5, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
  scale_x_continuous(breaks=c(1,20,40,60,80,100),labels=c(0,0.2,0.4,0.6,0.8,1),expand = c(0,0)) + 
  scale_y_continuous(breaks=c(1,20,40,60,80,100),labels=c(0,0.2,0.4,0.6,0.8,1),expand = c(0,0)) +
  coord_cartesian(xlim = c(1,100), ylim=c(1,100)) +
  xlab(label = "t (or s)") + ylab(label = "t' (or s')") +
  scale_fill_gradient(low='mediumblue', high ='yellow',
                      guide = guide_colorbar(barwidth=1, barheight = 18.5))
mine.heatmap # save as 5x4 inch


X <- Dat$B
X <- melt(X)
mine.heatmap <- ggplot(data = X, mapping = aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color="white") + theme_minimal() +
  theme(text=element_text(size=14, family="serif"), plot.title = element_text(hjust = 0.5, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
  scale_x_continuous(breaks=1:9,labels=0:8) +
  scale_y_continuous(breaks=1:9,labels=0:8) +
  coord_cartesian(xlim = c(1,9), ylim=c(1,9)) +
  xlab(label = "|t-t'| (discrete)") + ylab(label = "|s-s'| (discrete)") +
  scale_fill_gradient(low='mediumblue', high ='yellow',
                      guide = guide_colorbar(barwidth=1, barheight = 18.5))
mine.heatmap

Dat <- create_data(seed=17,buf=4,A1="brownian",A2="legendre",mask="ar",tau=0.5)
X <- Dat$B
X <- melt(X)
mine.heatmap <- ggplot(data = X, mapping = aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color="white") + theme_minimal() +
  theme(text=element_text(size=14, family="serif"), plot.title = element_text(hjust = 0.5, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
  scale_x_continuous(breaks=1:9,labels=0:8) +
  scale_y_continuous(breaks=1:9,labels=0:8) +
  coord_cartesian(xlim = c(1,9), ylim=c(1,9)) +
  xlab(label = "|t-t'| (discrete)") + ylab(label = "|s-s'| (discrete)") +
  scale_fill_gradient(low='mediumblue', high ='yellow',
                      guide = guide_colorbar(barwidth=1, barheight = 18.5))
mine.heatmap
