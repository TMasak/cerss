setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("library_smooth.R")
source("library_kronPCA.R")
setwd("./data")

process_data <- function(name,sim=25){
  load(paste(name,".RData",sep=""))
  ERR <- array(0,c(sim,7,7))
  for(i in 1:sim){
    ERR[i,,]<- Data[[i]]
  }
  ERR <- apply(ERR,c(2,3),mean)
  return(ERR)
}

process_cv <- function(name,SIM=25){
  load(paste(name,".RData",sep=""))
  ERR <- array(0,c(4,7))
  for(sim in 1:SIM){
    Akt <- Data[[sim]]
    for(pp in 1:7){
      ERR[1,pp] <- ERR[1,pp] + Akt[2,pp]
      ERR[2,pp] <- ERR[2,pp] + Akt[3,pp]
      iter <- which.min(Akt[6:10,pp])
      ERR[3,pp] <- ERR[3,pp] + Akt[iter,pp]
      ERR[4,pp] <- ERR[4,pp] + min(Akt[which(Akt[1:5,pp]!=0),pp])
      # there can be zeros due to convergence and algorithm stopping early
    }
  }
  ERR <- ERR/SIM
  return(ERR)
}

mycol <- c("blue","forestgreen","red","saddlebrown")
mypch <- c(0,1,4,3)
mylty <- rep(1,6) #c(4,2,3,1,5,6)
mylegend <- c("4D smoothing","one-step","proposed","BSA")
mylwd <- 2
mycex <- 1.2

library(plot3D)

##############
### Errors ###
##############

### Fourier

ERR <- process_data("sparse_fourier_data",100)
op <- par(mfrow=c(1,1),mar = c(3.2, 3, 0.2, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(2:7,ERR[1,2:7],type="b",ylim=c(0,1),xaxt="n",xlab="percentage observed",
     ylab="relative error",col="white")
axis(1,at=2:7,labels = c(2,5,10,20,40,70))
abline(h=seq(0,1,by=0.2),col="gray60",lty=5)
abline(h=ERR[7,7],lwd=mylwd,cex=mycex)
points(2:7,ERR[1,2:7],type="b",col=mycol[1],pch=mypch[1],lwd=mylwd,cex=mycex)
points(2:7,ERR[2,2:7],type="b",col=mycol[2],pch=mypch[2],lwd=mylwd,cex=mycex)
# points(2:4,EMP[4,2:4],type="b",col=mycol[3],pch=mypch[3],lwd=mylwd,cex=mycex)
points(2:4,ERR[4,2:4],type="b",col=mycol[3],pch=mypch[3],lwd=mylwd,cex=mycex)
legend("topright",legend=mylegend,col=c(mycol[3],mycol[1],mycol[2],"black"),pch=c(mypch[c(3,1,2)],NA),lty=c(NA,NA,NA,1),
       lwd=mylwd,y.intersp = 1.3,pt.cex=mycex) # save as 5x4


### Gneiting

ERR <- process_data("sparse_gneiting_data",100)
op <- par(mfrow=c(1,1),mar = c(3.2, 3, 0.2, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(2:7,ERR[1,2:7],type="b",ylim=c(0,1),xaxt="n",xlab="percentage observed",
     ylab="relative error",col="white")
axis(1,at=2:7,labels = c(2,5,10,20,40,70))
abline(h=seq(0,1,by=0.2),col="gray60",lty=5)
abline(h=ERR[7,7],lwd=mylwd,cex=mycex)
points(2:7,ERR[1,2:7],type="b",col=mycol[1],pch=mypch[1],lwd=mylwd,cex=mycex)
points(2:7,ERR[2,2:7],type="b",col=mycol[2],pch=mypch[2],lwd=mylwd,cex=mycex)
points(2:4,ERR[4,2:4],type="b",col=mycol[3],pch=mypch[3],lwd=mylwd,cex=mycex)
legend("topright",legend=mylegend,col=c(mycol[3],mycol[1],mycol[2],"black"),pch=c(mypch[c(3,1,2)],NA),lty=c(NA,NA,NA,1),
       lwd=mylwd,y.intersp = 1.3,pt.cex=mycex) # save as 5x4

### Brownian

ERR <- process_data("sparse_brownian_data",100)
op <- par(mfrow=c(1,1),mar = c(3.2, 3, 0.2, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(2:7,ERR[1,2:7],type="b",ylim=c(0,1),xaxt="n",xlab="percentage observed",
     ylab="relative error",col="white")
axis(1,at=2:7,labels = c(2,5,10,20,40,70))
abline(h=seq(0,1,by=0.2),col="gray60",lty=5)
abline(h=ERR[7,7],lwd=mylwd,cex=mycex)
points(2:7,ERR[1,2:7],type="b",col=mycol[1],pch=mypch[1],lwd=mylwd,cex=mycex)
points(2:7,ERR[2,2:7],type="b",col=mycol[2],pch=mypch[2],lwd=mylwd,cex=mycex)
points(2:4,ERR[4,2:4],type="b",col=mycol[3],pch=mypch[3],lwd=mylwd,cex=mycex)
legend("topright",legend=mylegend,col=c(mycol[3],mycol[1],mycol[2],"black"),pch=c(mypch[c(3,1,2)],NA),lty=c(NA,NA,NA,1),
       lwd=mylwd,y.intersp = 1.3,pt.cex=mycex)

### Mix - Fourier Legendre

ERR <- process_data("mix_fourier_legendre_data",100)
op <- par(mfrow=c(1,1),mar = c(3.2, 3, 0.2, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(2:7,ERR[1,2:7],type="b",ylim=c(0,1),xaxt="n",xlab="percentage observed",
     ylab="relative error",col="white")
axis(1,at=2:7,labels = c(2,5,10,20,40,70))
abline(h=seq(0,1,by=0.2),col="gray60",lty=5)
abline(h=ERR[7,7],lwd=mylwd,cex=mycex)
points(2:7,ERR[1,2:7],type="b",col=mycol[1],pch=mypch[1],lwd=mylwd,cex=mycex)
points(2:7,ERR[2,2:7],type="b",col=mycol[2],pch=mypch[2],lwd=mylwd,cex=mycex)
points(2:4,ERR[4,2:4],type="b",col=mycol[3],pch=mypch[3],lwd=mylwd,cex=mycex)
legend("bottomleft",legend=mylegend,col=c(mycol[3],mycol[1],mycol[2],"black"),pch=c(mypch[c(3,1,2)],NA),lty=c(NA,NA,NA,1),
       lwd=mylwd,y.intersp = 1.3,pt.cex=mycex, ncol=2, text.width=c(1.8,1.8))

########################
### Runtimes Fourier ###
########################

ERR <- process_data("runtimes_our_data")
EMP <- process_data("runtimes_empirical_data")
NNP <- process_data("runtimes_nonpooled_data")
op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(2:7,ERR[1,2:7],type="b",ylim=c(0,log(2560)),xaxt="n",yaxt="n",xlab="percentage observed",
     ylab="runtime [sec]",col="white")
axis(1,at=2:7,labels = c(2,5,10,20,40,70))
axis(2,at=c(0, 1.945910, 3.912023, 5.991465, 8.006368),labels=c(1,7,50,400,3000))
abline(h=c(0, 1.945910, 3.912023, 5.991465, 8.006368),col="gray60",lty=5)
points(2:6,log(NNP[4,2:6]),type="b",col="darkorange",pch=2,lwd=mylwd,cex=mycex)
points(2:7,log(ERR[4,2:7]),type="b",col=mycol[2],pch=mypch[2],lwd=mylwd,cex=mycex)
points(2:4,log(EMP[7,2:4]),type="b",col=mycol[3],pch=mypch[3],lwd=mylwd,cex=mycex)
legend(4.96,log(300),legend=c("4D smoothing","non-pooled data","proposed"),col=c(mycol[3],"darkorange",mycol[2]),
       pch=c(mypch[3],2,mypch[2]), lwd=mylwd, lty=c(NA,NA,NA),y.intersp = 1.3,pt.cex=mycex)

EMP[7,4]/ERR[4,4] # speed-up of separable model compared to the empirical estimator at the edge of feasibility

####################
### 4D bandwidth ###
####################

ERR <- array(0,c(2,4))
akt <- process_data("sparse_fourier_data")
ERR[1,1] <- akt[4,2]
akt <- process_data("sparse_fourier_empirical_CV_data")
ERR[2,1] <- akt[4,2]
akt <- process_data("sparse_brownian_data")
ERR[1,2] <- akt[4,2]
akt <- process_data("sparse_brownian_empirical_CV_data")
ERR[2,2] <- akt[4,2]
akt <- process_data("sparse_gneiting_data")
ERR[1,3] <- akt[4,2]
akt <- process_data("sparse_gneiting_empirical_CV_data")
ERR[2,3] <- akt[4,2]
akt <- process_data("mix_fourier_legendre_data")
ERR[1,4] <- akt[4,2]
akt <- process_data("mix_fourier_legendre_empirical_CV_data")
ERR[2,4] <- akt[4,2]

ERR

##########################
### 3rd iteration + CV ###
##########################

ERR <- array(0,c(12,6))
akt <- process_cv("sparse_fourier_data_cv",100)
ERR[1,] <- akt[1,2:7]
ERR[2,] <- akt[2,2:7]
ERR[3,] <- akt[3,2:7]
akt <- process_cv("sparse_gneiting_data_cv",100)
ERR[4,] <- akt[1,2:7]
ERR[5,] <- akt[2,2:7]
ERR[6,] <- akt[3,2:7]
akt <- process_cv("sparse_brownian_data_cv",100)
ERR[7,] <- akt[1,2:7]
ERR[8,] <- akt[2,2:7]
ERR[9,] <- akt[3,2:7]
akt <- process_cv("mix_fourier_legendre_data_cv",100)
ERR[10,] <- akt[1,2:7]
ERR[11,] <- akt[2,2:7]
ERR[12,] <- akt[3,2:7]

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 0.2, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(2:7,ERR[1,],type="b",ylim=c(0,1.1),xaxt="n",xlab="percentage observed",
     ylab="relative error",col="white")
axis(1,at=2:7,labels = c(2,5,10,20,40,70))
abline(h=seq(0,1,by=0.2),col="gray60",lty=5)
points(2:7,ERR[1,],type="b",col=mycol[1],pch=mypch[2],lwd=1.5,cex=mycex)
points(2:7,ERR[2,],type="b",col=mycol[1],pch=mypch[4],lwd=1.5,cex=mycex, lty=5)
points(2:7,ERR[3,],type="b",col=mycol[1],pch=mypch[3],lwd=1.5,cex=mycex, lty=3)
points(2:7,ERR[4,],type="b",col=mycol[2],pch=mypch[2],lwd=1.5,cex=mycex)
points(2:7,ERR[5,],type="b",col=mycol[2],pch=mypch[4],lwd=1.5,cex=mycex, lty=5)
points(2:7,ERR[6,],type="b",col=mycol[2],pch=mypch[3],lwd=1.5,cex=mycex, lty=3)
points(2:7,ERR[7,],type="b",col=mycol[3],pch=mypch[2],lwd=1.5,cex=mycex)
points(2:7,ERR[8,],type="b",col=mycol[3],pch=mypch[4],lwd=1.5,cex=mycex, lty=5)
points(2:7,ERR[9,],type="b",col=mycol[3],pch=mypch[3],lwd=1.5,cex=mycex, lty=3)
points(2:7,ERR[10,],type="b",col=mycol[4],pch=mypch[2],lwd=1.5,cex=mycex)
points(2:7,ERR[11,],type="b",col=mycol[4],pch=mypch[4],lwd=1.5,cex=mycex, lty=5)
points(2:7,ERR[12,],type="b",col=mycol[4],pch=mypch[3],lwd=1.5,cex=mycex, lty=3)

legend(4.75,0.98,legend=c("Fourier-Legendre","Fourier","Gneiting","Brownian"),
       col=mycol[c(4,1,2,3)],pt.lwd=mylwd,pt.cex=mycex,y.intersp = 1.3,lwd=1.5)
legend("topright",legend=c("2nd step    ","3rd step","CV"),pch=c(mypch[2],mypch[4],mypch[3]),ncol=3,
       pt.lwd=1.5,pt.cex=mycex)

########################
### simulation setup ###
########################

N <- 100       # no. of curves
K1 <- 20; K2 <- 20 # grid sizes 
ERR <- array(0,c(7,7)) # iter=1,2,3, runtime(for 2), bandwidths , fully observed x different percentages observed
C1 <- get_cov(K1,"fourier")
C1 <- C1/frobenius(C1)
C2 <- get_cov(K1,"brownian")
C2 <- C2/frobenius(C2)
C3 <- get_cov(K1,"legendre")
C3 <- C3/frobenius(C3)
C <- get_gneiting_cov(K1,K2,0.7) # for the purposes of error calculation
C <- C/frobenius(C)

C_mat <- tensor2matrix(C)
C_mat_perm <- VanLoansPerm(C_mat,K1,K2)
SVD <- svd(C_mat_perm)
Chat1 <- matrix(sign(SVD$u[1,1])*SVD$u[,1],ncol=K1)
Chat2 <- matrix(sign(SVD$v[1,1])*SVD$v[,1],ncol=K1)

op <- par(mfrow=c(1,1),c(0.5,0.5,0.5,2.5), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
persp3D(z = Chat1, theta = -30, phi = 30, border="black", lwd=0.5, bty = "g", expand=0.5, xlab="t", ylab="t'", zlab="value", contour=T)
persp3D(z = Chat2, theta = -30, phi = 30, border="black", lwd=0.5, bty = "g", expand=0.5, xlab="t", ylab="t'", zlab="value", contour=T)

# persp3D(z = C1, theta = -30, phi = 30, border="black", lwd=0.5, bty = "g", expand=0.5, xlab="t", ylab="t'", zlab="value", contour=T)
# persp3D(z = C2, theta = -30, phi = 30, border="black", lwd=0.5, bty = "g", expand=0.5, xlab="t", ylab="t'", zlab="value", contour=T)
# persp3D(z = C3, theta = -30, phi = 30, border="black", lwd=0.5, bty = "g", expand=0.5, xlab="t", ylab="t'", zlab="value", contour=T)

# shared limits for colkey for all the plots
op <- par(mfrow=c(1,1),mar = c(0.5,0,0,0), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=21, family='serif')
persp3D(z = C1, theta = -30, phi = 30, border="black", lwd=0.5, bty = "g", expand=0.5, xlab="t", ylab="t'", zlab="value", contour=T,
        colkey=F, clim=c(-0.0765,0.18))
persp3D(z = C2, theta = -30, phi = 30, border="black", lwd=0.5, bty = "g", expand=0.5, xlab="t", ylab="t'", zlab="value", contour=T,
        colkey=F, clim=c(-0.0765,0.18))
persp3D(z = C3, theta = -30, phi = 30, border="black", lwd=0.5, bty = "g", expand=0.5, xlab="t", ylab="t'", zlab="value", contour=T,
        colkey=F, clim=c(-0.0765,0.18))

# legend - save as 12x5 inch and then crop the pdf
op <- par(mfrow=c(1,1),mar = c(0.5,0,0,0), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
colkey(clim=c(-0.0765,0.19), side=1)
