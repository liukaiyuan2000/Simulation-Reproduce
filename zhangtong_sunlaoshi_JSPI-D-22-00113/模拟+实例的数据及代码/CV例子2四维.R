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
# 内积
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
  # g_l的方差
  const=(c(1:50)-0.5)*pi
  g=mvrnorm(n,rep(0,50),diag(const^(-2)))
  # tm是离散时间点，seq(0.005, 0.995, length.out = 100)
  phi=sqrt(2)*sin(const%*%t(tm))##50 times 100 matrix
  x=g%*%phi## n times 100 matrix
  z1=rnorm(n,1,1)
  z2=runif(n,0,pi)
  z3=rnorm(n,1,1)
  z4=runif(n,0,pi)
  z=cbind(z1,z2,z3,z4)
  u=mvrnorm(n,rep(0,4),covmatrix)
  w=z+u
  y=rep(0,n)
  error=rnorm(n,0,sigma)## sd 
  for(i in c(1:n)){
    y[i]=2*g[i,1]-g[i,2]+z[i,]%*%beta+error[i]+delta*exp(abs(z4[i]-1))*2/3
  }
  return(list(Y=y,X=x,Z=z,W=w))
}

# m是保留的主成分数
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
{
ztvalue=0.4
{
  tgridnum=100
  tm=c(1:tgridnum)*0.01-0.005
  alpha=2*sqrt(2)*sin(pi*tm/2)-sqrt(2)*sin(3*pi*tm/2)
  beta=c(1,0.8,0.5,0.3)
  sigma=0.5
  n=100
  diaguu=diag(c(ztvalue,ztvalue,ztvalue,0))
  delta=0
  rep=1000
}
cvresult=matrix(rep(0,10*rep),rep,10)
for(i in c(1:rep)){
  for(m in c(1:10)){
    for(j in c(1:n)){
      data=generate(n,diaguu,sigma,delta)
      x=data$X[j,]
      y=data$Y[j]
      w=data$W[j,]
      # 删一版本
      datax=data$X[-j,]
      datay=data$Y[-j]
      dataw=data$W[-j,]
      H=pca(datax,m,tm)$h
      eigf=pca(datax,m,tm)$eig.functions
      S=H%*%solve(t(H)%*%H)%*%t(H)
      # (4)式
      betahat=solve(t(dataw)%*%(diag(rep(1,(n-1)))-S)%*%dataw-(n-1)*diaguu)%*%t(dataw)%*%(diag(rep(1,(n-1)))-S)%*%datay
      alphahat=solve(t(H)%*%H)%*%t(H)%*%(datay-dataw%*%betahat)
      alphahatf=eigf%*%alphahat
      cvresult[i,m]=cvresult[i,m]+(y-t(w)%*%betahat-inpro(alphahatf,x,tm))^2-t(betahat)%*%diaguu%*%betahat
    }
  }
  if(i%%50==0){
    print(i)
  }
}


#write.csv(cvresult,"E:/first/eg2-44-cv.xls")
cvresult2=read.csv("E:/first/eg2-88-cv.xls")[,-1]


plot(cvresult2[,1]/100)
order(cvresult2[,2]/100)

ztindex=c(180,87,390)

ztvalue=0.8
{
  tgridnum=100
  tm=c(1:tgridnum)*0.01-0.005
  alpha=2*sqrt(2)*sin(pi*tm/2)-sqrt(2)*sin(3*pi*tm/2)
  beta=c(1,0.8,0.5,0.3)
  sigma=0.5
  n=100
  diaguu=diag(c(ztvalue,ztvalue,ztvalue,0))
  delta=0
  rep=1000
}
cvresult2=matrix(rep(0,10*rep),rep,10)

cvresult2[ztindex,]=rep(0,10)

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
      cvresult2[i,m]=cvresult2[i,m]+(y-t(w)%*%betahat-inpro(alphahatf,x,tm))^2-t(betahat)%*%diaguu%*%betahat
    }
  }
  #if(i%%50==0){
    print(i)
  #}
}
}

#write.csv(cvresult2,"E:/first/eg2-88-cv.xls")
# aaa=read.csv("E:/first/eg1-88-cv.xls")



plot(cvresult2[,1]/100)
plot(cvresult2[,2]/100)
plot(cvresult2[,3]/100)
plot(cvresult2[,4]/100)
plot(cvresult2[,5]/100)
plot(cvresult2[,6]/100)
plot(cvresult2[,7]/100)
plot(cvresult2[,8]/100)
plot(cvresult2[,9]/100)
plot(cvresult2[,10]/100)
order(cvresult2[,9]/100)
(apply(cvresult2/100,2,mean))
(apply(cvresult2/100,2,sd))
