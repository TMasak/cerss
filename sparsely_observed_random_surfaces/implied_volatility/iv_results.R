setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

mycol <- c("blue","forestgreen","red","saddlebrown")
mypch <- c(0,1,4,3)
mylty <- rep(1,6) #c(4,2,3,1,5,6)
mylegend <- c("4D smoothing","one-step","proposed","BSA")
mylwd <- 2
mycex <- 1.2

source("library_smooth.R")
library(plot3D)

load("TRout.RData")

Dell <- melt(TRout$Dell)
Dell <- Dell[complete.cases(Dell),]
Qualcomm <- melt(TRout$Qualcomm)
Qualcomm <- Qualcomm[complete.cases(Qualcomm),]
rozmezi <- c(min(Qualcomm),max(TRout$surface_mean))

op <- par(mfrow=c(1,1),mar = c(0.5,0,0,0), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=21, family='serif')
persp3D(z = TRout$surface_mean, theta = 40, phi = 11.5, border="black", lwd=0.5, bty = "g", expand=0.75, xlab="s (log-moneyness)",
        ylab="t (time to expiration)", zlab="log of implied volatility", contour=T, zlim=rozmezi, clim=rozmezi, colkey=F)
 # text(-0.37,-0.3,"s")

scatter3D(x=Dell$Var1, y=Dell$Var2, z=Dell$value, theta = 40, phi = 11.5,border="black", bty = "g", expand=0.75, pch=19, cex=1.8,
          ylab="t (time to expiration)", xlab="s (log-moneyness)", zlab="log of implied volatility", type="h",
          xlim=c(0,50), ylim=c(0,50), zlim=rozmezi, clim=rozmezi, colkey=F)
scatter3D(x=Dell$Var1, y=Dell$Var2, z=rep(rozmezi[1],length(Dell$Var1)),col="black", add=T)
scatter3D(x=Qualcomm$Var1, y=Qualcomm$Var2, z=Qualcomm$value, theta = 40, phi = 11.5,border="black", bty = "g", expand=0.75, pch=19,
          cex=1.8, ylab="t (time to expiration)", xlab="s (log-moneyness)", zlab="log of implied volatility", type="h",
          xlim=c(0,50), ylim=c(0,50), zlim=rozmezi, clim=rozmezi, colkey=F)
scatter3D(x=Qualcomm$Var1, y=Qualcomm$Var2, z=rep(rozmezi[1],length(Qualcomm$Var1)),col="black", add=T)

# legend - save as 12x5 inch and then crop the pdf
op <- par(mfrow=c(1,1),mar = c(0.5,0,0,0), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
colkey(clim=rozmezi, side=1, width=1)

### Marginal Kernels and their eigendecompositions

A <- TRout$Res$A*sqrt(max(TRout$Res$B)/max(TRout$Res$A))
B <- TRout$Res$B/sqrt(max(TRout$Res$B)/max(TRout$Res$A))

op <- par(mfrow=c(1,1),mar = c(0.5,0.5,0.5,2.5), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
persp3D(z=A, theta = -30, phi = 30, border="black", lwd=0.5, bty = "g", expand=0.5,
        xlab="s", ylab="s'", zlab="value", contour=T)
persp3D(z=B, theta = 150, phi = 30, border="black", lwd=0.5, bty = "g", expand=0.5,
        xlab="t", ylab="t'", zlab="value", contour=T)

K1 <- K2 <- 50
op <- par(mfrow=c(1,1),mar = c(3.2, 2, 0, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')
plot(seq(0.5,1.5,length.out=K1),rep(0,K1),type="b",ylim=c(-0.4,0.5),xaxt="n",xlab="s (log-moneyness)",
     ylab="",col="white")
axis(1)
abline(h=seq(-0.4,0.4,by=0.2),col="gray60",lty=5)
points(seq(0.5,1.5,length.out=K1), TRout$EIG1$vectors[,1],col="red",type="l",lwd=mylwd)
points(seq(0.5,1.5,length.out=K1), TRout$EIG1$vectors[,2],col="blue",type="l",lwd=mylwd)
points(seq(0.5,1.5,length.out=K1), TRout$EIG1$vectors[,3],col="forestgreen",type="l",lwd=mylwd)
legend("bottomleft", ncol=3, legend= paste( round(TRout$EIG1$values[1:3]/sum(TRout$EIG1$values),2)*100 , "%"),
       col=c('red','blue','forestgreen'), lty=1 , lwd=mylwd, text.width=0.15*c(1,1,1))


plot(seq(14,365,length.out=K1),rep(0,K1),type="b",ylim=c(-0.4,0.5),xaxt="n",xlab="t (time to expiration)",
     ylab="",col="white")
axis(1)
abline(h=seq(-0.4,0.4,by=0.2),col="gray60",lty=5)
points(seq(14,365,length.out=K1), TRout$EIG2$vectors[,1],col="red",type="l",lwd=mylwd)
points(seq(14,365,length.out=K1), TRout$EIG2$vectors[,2],col="blue",type="l",lwd=mylwd)
points(seq(14,365,length.out=K1), TRout$EIG2$vectors[,3],col="forestgreen",type="l",lwd=mylwd)
legend("bottomleft", ncol=3, legend=c("96 %", "3 %", "0.5 %"),
       col=c('red','blue','forestgreen'), lty=1 , lwd=mylwd, text.width=35*c(1,1,1))

### Prediction plot

Dell_fit <- melt(TRout$Dell)
SURF <- melt(TRout$posterior$Mean)
Dell_fit <- SURF$value[complete.cases(Dell_fit)]

# scatter3D(x=Dell$Var1, y=Dell$Var2, z=Dell$value, theta = 40, phi = 11.5,border="black", bty = "b2", expand=0.75, pch=19, cex=1.8,
#           xlab="t", ylab="t'", zlab="value", xlim=c(0,50), ylim=c(0,50), zlim=rozmezi, clim=rozmezi,
#           surf=list(x=1:K1, y=1:K2, z=TRout$posterior$Mean, fit=Dell_fit, alpha=0.5))
# persp3D(x=1:K1, y=1:K2, z=TRout$posterior$Upper, add=T, col="black", alpha=0.1)
# persp3D(x=1:K1, y=1:K2, z=TRout$posterior$Lower, add=T, col="black", alpha=0.1)

op <- par(mfrow=c(1,1),mar = c(0.5,1.5,0.5,2), bty = "n",mgp=c(1.8,0.7,0), las=0,
          ps=14, family='serif')

Theta=-40
scatter3D(x=Dell$Var1, y=Dell$Var2, z=Dell$value, theta = Theta, phi = 11.5,border="black", bty = "b2", expand=0.75, pch=1, cex=1.7,
          xlab="s (log-moneyness)", ylab="t (time to expiration)", zlab="log of implied volatility", xlim=c(0,50), ylim=c(0,50),
          zlim=rozmezi, clim=rozmezi,surf=list(x=1:K1, y=1:K2, z=TRout$posterior$Mean, fit=Dell_fit, alpha=0), colkey=F, col="deeppink")
scatter3D(x=Dell$Var1, y=Dell$Var2, z=Dell$value, pch=19, cex=1.6, add=T, clim=rozmezi,colkey=F)
# text(-0.36,-0.3,"t")
subMean <- TRout$posterior$Mean[-seq(2,K1-2,by=2),-seq(2,K2-2,by=2)]
ribbon3D(x=c(seq(1,K1,by=2),K1) ,y=c(seq(1,K2,by=2),K2), z=subMean, along="xy", add=T, colkey=F, lwd=1, space=0.7)
persp3D(x=1:K1, y=1:K2, z=TRout$posterior$Upper, add=T, col="black", alpha=0.1)
persp3D(x=1:K1, y=1:K2, z=TRout$posterior$Lower, add=T, col="black", alpha=0.1)

Theta=40
scatter3D(x=Dell$Var1, y=Dell$Var2, z=Dell$value, theta = Theta, phi = 11.5,border="black", bty = "b2", expand=0.75, pch=1, cex=1.7,
          xlab="s (log-moneyness)", ylab="t (time to expiration)", zlab="log of implied volatility", xlim=c(0,50), ylim=c(0,50),
          zlim=rozmezi, clim=rozmezi, surf=list(x=1:K1, y=1:K2, z=TRout$posterior$Mean, fit=Dell_fit, alpha=0), colkey=F, col="deeppink")
scatter3D(x=Dell$Var1, y=Dell$Var2, z=Dell$value, pch=19, cex=1.6, add=T, clim=rozmezi,colkey=F)
subMean <- TRout$posterior$Mean[-seq(2,K1-2,by=2),-seq(2,K2-2,by=2)]
ribbon3D(x=c(seq(1,K1,by=2),K1) ,y=c(seq(1,K2,by=2),K2), z=subMean, along="xy", add=T, colkey=F, lwd=1, space=0.7)
persp3D(x=1:K1, y=1:K2, z=TRout$posterior$Upper, add=T, col="black", alpha=0.1)
persp3D(x=1:K1, y=1:K2, z=TRout$posterior$Lower, add=T, col="black", alpha=0.1)


# maybe display colorkey exponentiated?

### Boxplots

ratio <- TRout$ratio
estimator <- c(rep("a_separable",196),rep("a_4D",196),
               rep("b_separable",196),rep("b_4D",196),
               rep("c_separable",196),rep("c_4D",196),
               rep("d_separable",196),rep("d_4D",196),
               rep("e_separable",196),rep("e_4D",196))
Data <- data.frame(ratio=ratio,estimator=estimator)
Data <- Data[complete.cases(Data),]

op <- par(mfrow=c(1,1),mar = c(3.2, 3, 1.6, 0.2), bty = "n",mgp=c(1.8,0.7,0), las=0, ps=14, family='serif')
plot(1,1,col="white",xlim=c(1,14),ylim=c(-3.5,3.5),xaxt="n",yaxt="n",xlab="hold-out pattern",ylab="RMSE")
# axis(1,at=c(1,2),labels=c("(a)","(a)"))
axis(1,at=c(1,2),labels=c("",""))
axis(1,at=c(1.5),labels=c("(a)"),tick=F)
axis(1,at=c(4,5),labels=c("",""))
axis(1,at=c(4.5),labels=c("(b)"),tick=F)
axis(1,at=c(7,8),labels=c("",""))
axis(1,at=c(7.5),labels=c("(c)"),tick=F)
axis(1,at=c(10,11),labels=c("",""))
axis(1,at=c(10.5),labels=c("(d)"),tick=F)
axis(1,at=c(13,14),labels=c("",""))
axis(1,at=c(13.5),labels=c("(e)"),tick=F)
where <- c(1/27,1/9,1/3,1,3,9,27)
axis(2,at=log(where),labels=c("1/27","1/9","1/3","1","3","9","27"),lwd.ticks=c(1,1,1.5,1,1),las=1)
axis(2,at=0,labels=1,lwd=1.5,las=1)
abline(h=log(where),col="gray60")
abline(h=0,lwd=1.5)
boxplot(log(ratio) ~ estimator, data=Data, at=c(1,2,4,5,7,8,10,11,13,14),add=T,xaxt="n",yaxt="n",
        col=c("gray55","gray90"))#,border=c("gray40","gray80"))

legend("topright", legend = c("4D smoothing","proposed"),pch=c(22,22),pt.cex = 2,pt.bg=c("gray50","gray90"))

median_sep <- c(median(Data$ratio[Data$estimator=="a_separable"]),
                median(Data$ratio[Data$estimator=="b_separable"]),
                median(Data$ratio[Data$estimator=="c_separable"]),
                median(Data$ratio[Data$estimator=="d_separable"]),
                median(Data$ratio[Data$estimator=="e_separable"]))
median_4D <- c(median(Data$ratio[Data$estimator=="a_4D"]),
               median(Data$ratio[Data$estimator=="b_4D"]),
               median(Data$ratio[Data$estimator=="c_4D"]),
               median(Data$ratio[Data$estimator=="d_4D"]),
               median(Data$ratio[Data$estimator=="e_4D"]))
text(c(1,4,7,10,13)+1,log(median_sep)+0.26,labels=round(median_sep,digit=2))
text(c(1,4,7,10,13),log(median_4D)+0.26,labels=round(median_4D,digit=2),col="white")
