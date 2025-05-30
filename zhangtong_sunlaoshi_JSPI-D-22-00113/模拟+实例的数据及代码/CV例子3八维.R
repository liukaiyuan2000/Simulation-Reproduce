install.packages("MASS")
install.packages("fda.usc")
install.packages("latex2exp")
library(MASS)
library(fda.usc)
library(latex2exp)
#kernel function
fkernel=function(x,h){
  return(dnorm(abs(x)/h))
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
#generate data function
generate=function(n,covmatrix,sigma,delta){
  const=(c(1:50)-0.5)*pi
  g=mvrnorm(n,rep(0,50),diag(const^(-2)))
  phi=sqrt(2)*sin(const%*%t(tm))##50 times 100 matrix
  x=g%*%phi## n times 100 matrix
  z1=rnorm(n,1,1)
  z2=rnorm(n,1,1)
  z3=rnorm(n,1,1)
  z4=rnorm(n,1,1)
  z5=runif(n,0,pi)
  z6=runif(n,0,pi)
  z7=runif(n,0,pi)
  z8=runif(n,0,pi)
  z=cbind(z1,z2,z3,z4,z5,z6,z7,z8)
  u=mvrnorm(n,rep(0,8),covmatrix)
  w=z+u
  y=rep(0,n)
  error=rnorm(n,0,sigma)## sd 
  for(i in c(1:n)){
    y[i]=g[i,1]+g[i,2]+z[i,]%*%beta+error[i]#+delta/8*z[i,3]^3
  }
  return(list(Y=y,X=x,Z=z,W=w))
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

##CV PROCODURE

ztvalue=0.8

  {
    tgridnum=100
    tm=c(1:tgridnum)*0.01-0.005
    #alpha=sqrt(2)*sin(pi*tm/2)+sqrt(2)*sin(3*pi*tm/2)
    beta=c(0.8,0.6,0.4,0.2,0.8,0.6,0.4,0.2)
    sigma=0.5
    n=100
    diaguu=diag(c(ztvalue,ztvalue,ztvalue,ztvalue,ztvalue,ztvalue,0,0))
    delta=0
    rep=1000
  }
  cvresult=matrix(rep(0,10*rep),rep,10)
  
  
  ztindex=c(239,586,943,83,370)
  cvresult[ztindex,]=rep(0,10)
  
  for(i in ztindex){
    for(m in c(1:10)){
      for(j in c(1:n)){
        data=generate(n,diaguu,sigma,delta)
        x=data$X[j,]
        y=data$Y[j]
        w=data$W[j,]
        datax=data$X[-j,]
        datay=data$Y[-j]
        dataw=data$W[-j,]
        H=pca(datax,m,tm)$h
        eigf=pca(datax,m,tm)$eig.functions
        S=H%*%solve(t(H)%*%H)%*%t(H)
        betahat=solve(t(dataw)%*%(diag(rep(1,(n-1)))-S)%*%dataw-(n-1)*diaguu)%*%t(dataw)%*%(diag(rep(1,(n-1)))-S)%*%datay
        alphahat=solve(t(H)%*%H)%*%t(H)%*%(datay-dataw%*%betahat)
        alphahatf=eigf%*%alphahat
        cvresult[i,m]=cvresult[i,m]+(y-t(w)%*%betahat-inpro(alphahatf,x,tm))^2-t(betahat)%*%diaguu%*%betahat
      }
    }
    #if(i%%100==0){
      print(i)
    #}
  }
  print(apply(cvresult,2,mean)/100)
  plot(apply(cvresult,2,mean)/100)
  
  max(cvresult)
boxplot(cvresult)
cvresult[1,]
quantile(cvresult,0.99)

mean(cvresult[-71,1])
plot(cvresult[-71,1])

mean(cvresult[intersect(which(cvresult[,1]<=quantile(cvresult[,1],0.99)),which(cvresult[,1]>=quantile(cvresult[,1],0.01))),1])/100
mean(cvresult[intersect(which(cvresult[,2]<=quantile(cvresult[,2],0.99)),which(cvresult[,2]>=quantile(cvresult[,2],0.01))),2])/100
mean(cvresult[intersect(which(cvresult[,3]<=quantile(cvresult[,3],0.99)),which(cvresult[,3]>=quantile(cvresult[,3],0.01))),3])/100
mean(cvresult[intersect(which(cvresult[,4]<=quantile(cvresult[,4],0.99)),which(cvresult[,4]>=quantile(cvresult[,4],0.01))),4])/100
mean(cvresult[intersect(which(cvresult[,5]<=quantile(cvresult[,5],0.99)),which(cvresult[,5]>=quantile(cvresult[,5],0.01))),5])/100
mean(cvresult[intersect(which(cvresult[,6]<=quantile(cvresult[,6],0.99)),which(cvresult[,6]>=quantile(cvresult[,6],0.01))),6])/100
mean(cvresult[intersect(which(cvresult[,7]<=quantile(cvresult[,7],0.99)),which(cvresult[,7]>=quantile(cvresult[,7],0.01))),7])/100
mean(cvresult[intersect(which(cvresult[,8]<=quantile(cvresult[,8],0.99)),which(cvresult[,8]>=quantile(cvresult[,8],0.01))),8])/100
mean(cvresult[intersect(which(cvresult[,9]<=quantile(cvresult[,9],0.99)),which(cvresult[,9]>=quantile(cvresult[,9],0.01))),9])/100
mean(cvresult[intersect(which(cvresult[,10]<=quantile(cvresult[,10],0.99)),which(cvresult[,10]>=quantile(cvresult[,10],0.01))),10])/100



#write.csv(cvresult,"E:/first/eg3-88-cv-new.xls")


plot(cvresult[,1])
order(cvresult[,1])

plot(apply(cvresult/100, 2,mean))

#cvresult=read.csv("E:/first/eg3-88-cv-new.xls")[,-1]
#str(aaa)
aaa=aaa[,-1]
#str(aaa)
print(apply(aaa,2,mean)/100)
resaaa=c()
for (i in c(1:10)) {
  resaaa[i]=mean(aaa[intersect(which(aaa[,i]<=quantile(aaa[,i],0.9)),which(aaa[,i]>=quantile(aaa[,i],0.1))),i])/100
}
plot(resaaa)
print(resaaa)


plot(cvresult[,1])
plot(sort(cvresult[,1]))


for (i in c(1:10)) {
   print(mean(cvresult[which(abs(cvresult[,i])<1000),i])/100)
}
diaguu
for (i in c(1:10)) {
  print(mean(cvresult[intersect(which(cvresult[,i]<=1050),which(cvresult[,i]>=-850)),i])/100)
}
mean(cvresult[which(abs(cvresult[,1])<1000),1])/100
mean(cvresult[which(abs(cvresult[,2])<1000),2])/100
mean(cvresult[which(abs(cvresult[,3])<1000),3])/100


plot(cvresult[which(abs(cvresult[,1])<1000),1])
order(cvresult[,1])

cvresult[742,1]

plot(aaa[,10])

