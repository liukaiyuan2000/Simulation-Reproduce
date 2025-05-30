install.packages("MASS")
install.packages("fda.usc")
install.packages("latex2exp")
library(MASS)
library(fda.usc)
library(latex2exp)
#real data analysis part
#read data and pretreatment
doubleindex1=c(12,13,15,16,17,26,28,36,37,38,39,40,58,59,61,64,75,86,89,139,145,186)
doubleindex2=c(48,188,51,54,190,138,29,110,116,117,142,123,135,192,62,150,136,180,181,204,176,215)
deleterestindex=setdiff(c(1:215),union(c(34,35,43,44,45,129,140,172,186,185,99,139),doubleindex2))
nnew=length(deleterestindex)

{
  data=read.table("D://tecatornew.txt")
  xreal=matrix(,nnew,100)
  yreal=c()
  zreal=matrix(,nnew,2)
  ii=1
  for(i in deleterestindex){
    for(j in c(1:20)){
      xreal[ii,((5*j-4):(5*j))]=as.numeric(data[(25*i-25+j),])
    }
    yreal[ii]=as.numeric(data[25*i,4])#fat
    zreal[ii,1]=as.numeric(data[25*i,5])#moisture
    zreal[ii,2]=as.numeric(data[25*i,3])#pr
    ii=ii+1
  }
tmreal=c(1:100)*2+849
}

#centerlized
{
x1c=matrix(,nnew,100)
 for(i in c(1:100)){
  x1c[,i]=xreal[,i]-mean(xreal[,i])
 }
yreal=yreal-mean(yreal)
 for(i in c(1:2)){
  zreal[,i]=zreal[,i]-mean(zreal[,i])
 }
}
#plot the functional data after pretreatment
par(mfrow=c(1,2))
plot(tmreal,xreal[1,],type = "l",ylim = c(2,5.5),col="grey",xlab= "Wavelength ",ylab = "Absorbance")
for (i in c(1:nnew)){
  par(new=TRUE)
  plot(tmreal,xreal[i,],type = "l",xlab= " ",ylab = "",ylim = c(2,5.5),col="grey")
}
#PCA function
pca=function(X,m,t){
  if(m>1){
    eig.functions=eigen(cov(X))$vector[,1:m]*sqrt(length(t))
    eig.scores=matrix(,dim(X)[1],m)
    for(i in c(1:dim(X)[1])){
      for(j in c(1:m)){
        eig.scores[i,j]=inpro(X[i,],eig.functions[,j],t)
      }
    }
    eig.values <- apply((eig.scores)^2,2,mean)
    
  }
  else if(m==1){
    eig.functions=eigen(cov(X))$vector[,1]*sqrt(length(t))
    eig.scores=rep(0,dim(X)[1])
    for(i in c(1:dim(X)[1])){
      eig.scores[i]=inpro(X[i,],eig.functions,t)
    }
    eig.values <- mean((eig.scores)^2)
  }
  return(list("eig.values"=eig.values,"eig.functions"=eig.functions,"h"=eig.scores))
}
#inpro function
inpro=function(x,y,t){
  x=as.numeric(x)
  y=as.numeric(y)
  n=length(x)
  delta=(t[n]-t[1])/(n-1)
  result=x%*%y*delta
  return(result)
}
#kernel function
fkernel=function(x,h){
  return(dnorm(abs(x)/h))
}


{
#参数选取
diaguu=diag(c(3,0))
mcv=3
#估计部分一些统计量的计�?
{ 
  H=pca(x1c,mcv,tmreal)$h
  eigf=pca(x1c,mcv,tmreal)$eig.functions
  S=H%*%solve(t(H)%*%H)%*%t(H)
  d=matrix(rep(0,nnew^2),nnew,nnew)
  for(i in c(2:nnew)){
    for(j in c(1:(i-1))){
      d[i,j]=sqrt(inpro(x1c[i,]-x1c[j,],x1c[i,]-x1c[j,],tmreal))
    }
  }
  wreal=zreal+mvrnorm(nnew,rep(0,2),diaguu)
  d=d+t(d)
  h1=sd(d)*1.06*nnew^(-0.2)
  hw22=sd(wreal[,2])*1.06*nnew^(-0.2)
  Kw=matrix(rep(0,nnew^2),nnew,nnew)
  for(i in c(2:nnew)){
    for(j in c(1:(i-1))){
      Kw[i,j]=fkernel(wreal[i,2]-wreal[j,2],hw22)*fkernel(d[i,j],h1)
    }
  }
  Kw=Kw+t(Kw)
}


u1=rnorm(nnew,0,sqrt(diaguu[1,1]))
u2=rnorm(nnew,0,sqrt(diaguu[1,1]))
diaguuhat=diag(c(mean((u1-u2)^2)/2,0))


###### modify diaguuhat to diaguu or 0, lead to known method or naive method
#estimation
{
  betahat=solve(t(wreal)%*%(diag(rep(1,nnew))-S)%*%wreal-nnew*diaguuhat)%*%t(wreal)%*%(diag(rep(1,nnew))-S)%*%yreal
  alphahat=solve(t(H)%*%H)%*%t(H)%*%(yreal-wreal%*%betahat)
  alphahatf=eigf%*%alphahat
  ehat=rep(0,nnew)
  for(i in c(1:nnew)){
    ehat[i]=yreal[i]-t(wreal[i,])%*%betahat-inpro(alphahatf,x1c[i,],tmreal)
  }
  statvalue=t(ehat)%*%Kw%*%ehat
}
##bootstrap
{
  ehatcun=ehat
  B=1000
  bvalue=c()
  for(index in c(1:B)){
    eb=c()
    yb=c()
    for(i in c(1:nnew)){
      if(runif(1)<(5+sqrt(5))/10){
        eb[i]=(1-sqrt(5))/2
      }
      else
        eb[i]=(1+sqrt(5))/2
    }
    ebnew=eb*ehatcun
    for(i in c(1:nnew)){
      yb[i]=t(wreal[i,])%*%betahat+inpro(alphahatf,x1c[i,],tmreal)+ebnew[i]
    }
    bbetahat=solve(t(wreal)%*%(diag(rep(1,nnew))-S)%*%wreal)%*%t(wreal)%*%(diag(rep(1,nnew))-S)%*%yb
    balphahat=solve(t(H)%*%H)%*%t(H)%*%(yb-wreal%*%bbetahat)
    balphahatf=eigf%*%balphahat
    ehat=rep(0,nnew)
    for(i in c(1:nnew)){
      ehat[i]=yb[i]-t(wreal[i,])%*%bbetahat-inpro(balphahatf,x1c[i,],tmreal)
    }
    bvalue[index]=t(ehat)%*%Kw%*%ehat
  }
} 
mean(rep(abs(statvalue),B)<abs(bvalue))

#####write.csv(cbind(x2,y2),"E:/firrda/30-unknown-221.csv")

{
x2=yreal-ehatcun
y2=ehatcun
x=x2
y=y2
h=density(x)$bw
#functions to plot 95\% CI band
fx.hat <- function(z, h) {
  dnorm(abs(z - x)/h)/h
}
NWSMOOTH <- function(h, y, x) {
  n <- length(y)
  s.hat <- rep(0, n)
  for (i in 1:n) {
    a <- fx.hat(x[i], h)
    s.hat[i] <- sum(y * a/sum(a))
  }
  return(s.hat)
}
NWsmooth.val <- NWSMOOTH(h, y, x)

fSMOOTH <- function(h, x) {
  n <- length(x)
  s.hat <- rep(0, n)
  for (i in 1:n) {
    a <- fx.hat(x[i], h)
    s.hat[i] <- mean(a)
  }
  return(s.hat)
}
fsmooth.val <- fSMOOTH(h, x)
#standand error
sdes=1.96*sqrt(var(y)*0.28/fsmooth.val/nnew/h)
#plot the figure 
xlab=expression(W^{T}*hat(beta)+integral(hat(alpha)(t)*X(t)*dt,0,1))
ylab=expression(hat(e))
plot(x, y,xlim=c(-20,30),ylim=c(-6,6),xlab=xlab,ylab=ylab,main = "(a)",col="blue",pch=1,cex=1.5)
abline(h=0,lty=2,lwd=3,col="black")
par(new=TRUE)
plot(x[order(x)],NWsmooth.val[order(x)],lty =1,lwd=2,col = 2,xlim=c(-20,30),ylim=c(-6,6),type = "l",xlab="",ylab="")
par(new=TRUE)
plot(x[order(x)],(NWsmooth.val+sdes)[order(x)],lty=2,lwd=2,col=2,xlim=c(-20,30),ylim=c(-6,6),type = "l",xlab="",ylab="")
par(new=TRUE)
plot(x[order(x)],(NWsmooth.val-sdes)[order(x)],lty=2,lwd=2,col=2,xlim=c(-20,30),ylim=c(-6,6),type = "l",xlab="",ylab="")

}




###plot the final figure
par(mfrow=c(1,2),mar=c(5,4,3,2)+0.1)


{
  {
###known method
x=read.csv("E:/firrda/30-known.csv")[,2]
y=read.csv("E:/firrda/30-known.csv")[,3]
h=density(x)$bw

fx.hat <- function(z, h) {
  dnorm(abs(z - x)/h)/h
}
NWSMOOTH <- function(h, y, x) {
  n <- length(y)
  s.hat <- rep(0, n)
  for (i in 1:n) {
    a <- fx.hat(x[i], h)
    s.hat[i] <- sum(y * a/sum(a))
  }
  return(s.hat)
}
NWsmooth.val <- NWSMOOTH(h, y, x)

fSMOOTH <- function(h, x) {
  n <- length(x)
  s.hat <- rep(0, n)
  for (i in 1:n) {
    a <- fx.hat(x[i], h)
    s.hat[i] <- mean(a)
  }
  return(s.hat)
}
fsmooth.val <- fSMOOTH(h, x)

sdes=1.96*sqrt(var(y)*0.28/fsmooth.val/nnew/h)


xlab=TeX('$\\mathbf{W}^{T} \\hat{\\beta}_n +\\int_{850}^{1050} \\hat{\\alpha}_n(t)X(t)dt$')
ylab=TeX('$\\hat{\\epsilon}_n$')
plot(x, y,xlim=c(-20,30),ylim=c(-6,6),xlab="",ylab="",main = "(a)",col="blue",pch=16,cex=1.2)
mtext(ylab, side = 2, line = 2,cex=1)
mtext(xlab, side = 1, line = 4,cex=1)
abline(h=0,lty=2,lwd=3,col="black")
par(new=TRUE)
plot(x[order(x)],NWsmooth.val[order(x)],lty =1,lwd=2,col = 2,xlim=c(-20,30),ylim=c(-6,6),type = "l",xlab="",ylab="")
par(new=TRUE)
plot(x[order(x)],(NWsmooth.val+sdes)[order(x)],lty=2,lwd=2,col=2,xlim=c(-20,30),ylim=c(-6,6),type = "l",xlab="",ylab="")
par(new=TRUE)
plot(x[order(x)],(NWsmooth.val-sdes)[order(x)],lty=2,lwd=2,col=2,xlim=c(-20,30),ylim=c(-6,6),type = "l",xlab="",ylab="")
}

####unknown method
  {
x=read.csv("E:/firrda/30-unknown-221.csv")[,2]
y=read.csv("E:/firrda/30-unknown-221.csv")[,3]

h=density(x)$bw

fx.hat <- function(z, h) {
  dnorm(abs(z - x)/h)/h
}
NWSMOOTH <- function(h, y, x) {
  n <- length(y)
  s.hat <- rep(0, n)
  for (i in 1:n) {
    a <- fx.hat(x[i], h)
    s.hat[i] <- sum(y * a/sum(a))
  }
  return(s.hat)
}
NWsmooth.val <- NWSMOOTH(h, y, x)

fSMOOTH <- function(h, x) {
  n <- length(x)
  s.hat <- rep(0, n)
  for (i in 1:n) {
    a <- fx.hat(x[i], h)
    s.hat[i] <- mean(a)
  }
  return(s.hat)
}
fsmooth.val <- fSMOOTH(h, x)

sdes=1.96*sqrt(var(y)*0.28/fsmooth.val/nnew/h)

xlab=TeX('$\\mathbf{W}^{T} \\hat{\\beta}_n +\\int_{850}^{1050} \\hat{\\alpha}_n(t)X(t)dt$')
ylab=TeX('$\\hat{\\epsilon}_n$')
plot(x, y,xlim=c(-20,30),ylim=c(-6,6),xlab="",ylab="",main = "(b)",col="blue",pch=16,cex=1.2)
mtext(ylab, side = 2, line = 2,cex=1)
mtext(xlab, side = 1, line = 4,cex=1)
abline(h=0,lty=2,lwd=3,col="black")
par(new=TRUE)
plot(x[order(x)],NWsmooth.val[order(x)],lty =1,lwd=2,col = 2,xlim=c(-20,30),ylim=c(-6,6),type = "l",xlab="",ylab="")
par(new=TRUE)
plot(x[order(x)],(NWsmooth.val+sdes)[order(x)],lty=2,lwd=2,col=2,xlim=c(-20,30),ylim=c(-6,6),type = "l",xlab="",ylab="")
par(new=TRUE)
plot(x[order(x)],(NWsmooth.val-sdes)[order(x)],lty=2,lwd=2,col=2,xlim=c(-20,30),ylim=c(-6,6),type = "l",xlab="",ylab="")
}

###naive method
  {
x=read.csv("E:/firrda/30-naive-026.csv")[,2]
y=read.csv("E:/firrda/30-naive-026.csv")[,3]
h=density(x)$bw

fx.hat <- function(z, h) {
  dnorm(abs(z - x)/h)/h
}
NWSMOOTH <- function(h, y, x) {
  n <- length(y)
  s.hat <- rep(0, n)
  for (i in 1:n) {
    a <- fx.hat(x[i], h)
    s.hat[i] <- sum(y * a/sum(a))
  }
  return(s.hat)
}
NWsmooth.val <- NWSMOOTH(h, y, x)

fSMOOTH <- function(h, x) {
  n <- length(x)
  s.hat <- rep(0, n)
  for (i in 1:n) {
    a <- fx.hat(x[i], h)
    s.hat[i] <- mean(a)
  }
  return(s.hat)
}
fsmooth.val <- fSMOOTH(h, x)

sdes=1.96*sqrt(var(y)*0.28/fsmooth.val/nnew/h)

xlab=TeX('$\\mathbf{W}^{T} \\hat{\\beta}^{naive}_n +\\int_0^1 \\hat{\\alpha}_n^{naive}(t)X(t)dt$')
ylab=TeX('$\\hat{\\epsilon}^{naive}_n$')
plot(x, y,xlim=c(-20,30),ylim=c(-6,6),xlab="",ylab="",main = "(c)",col="blue",pch=1,cex=1.5)
mtext(ylab, side = 2, line = 2,cex=1)
mtext(xlab, side = 1, line = 4,cex=1)
abline(h=0,lty=2,lwd=3,col="black")
par(new=TRUE)
plot(x[order(x)],NWsmooth.val[order(x)],lty =1,lwd=2,col = 2,xlim=c(-20,30),ylim=c(-6,6),type = "l",xlab="",ylab="")
par(new=TRUE)
plot(x[order(x)],(NWsmooth.val+sdes)[order(x)],lty=2,lwd=2,col=2,xlim=c(-20,30),ylim=c(-6,6),type = "l",xlab="",ylab="")
par(new=TRUE)
plot(x[order(x)],(NWsmooth.val-sdes)[order(x)],lty=2,lwd=2,col=2,xlim=c(-20,30),ylim=c(-6,6),type = "l",xlab="",ylab="")
}

}





