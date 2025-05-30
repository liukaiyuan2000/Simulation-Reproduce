library(survival)
n=50 
X=rnorm(n,2,1)
Ystar=2*X+3+rnorm(n)
Y=pmax(Ystar,quantile(Ystar,0.3))
C=min(Y)
delta=1*(Y==Ystar)
dat<-data.frame(X=X,Y=Y,delta=delta,Ystar=Ystar)  

pch_site<-5*(1-dat$delta)
pch_col<-2*(1-dat$delta)+2
plot(Ystar~X,dat,col=pch_col,pch=pch_site,xlim=c(-1,3.5),ylim=c(2,13),xlab="X",ylab="The observations")
abline(a=3,b=2,col=1,lwd=3) # the true line
fit=lm(Ystar~X)
abline(fit,col=2,lwd=3)  #the least square estimator
par(new=TRUE)
pch_site<-20*dat$delta+1
pch_col<-2*(1-dat$delta)+3
plot(Y~X,dat,col=2,lty=2,pch=pch_site,xlim=c(-1,3.5),ylim=c(2,13),xlab="X",ylab="The observations")
fit2=lm(Y~X)
abline(fit2,col="deeppink",lty=6,lwd=4) # the naive least square estimator

fit3=survreg(Surv(Y, Y>C, type='left') ~X, dist='gaussian')
abline(fit3,col=4,lty=4,lwd=4)  #MLE 
Y1=Y[Y==Ystar]
X1=X[Y==Ystar]
fit4=lm(Y1~X1)
abline(fit4,col="blueviolet",lty=5,lwd=4) #2-LSE  
title("n=50", col.main = "blue")

legend(-0.5,13, c("Golden line","LSE(Ystar)","NLSE(Y)","MLE","2-LSE(delta=1)"), cex=0.8,lwd=3,col=c("black","red","deeppink","blue","blueviolet"),lty=c(1,1,6,4,5))

