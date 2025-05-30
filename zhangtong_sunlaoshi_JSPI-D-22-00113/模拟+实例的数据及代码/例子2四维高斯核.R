install.packages("MASS")
install.packages("fda.usc")
library(MASS)
library(fda.usc)
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
  z2=runif(n,0,pi)
  z3=rnorm(n,1,1)
  z4=runif(n,0,pi)
  z=cbind(z1,z2,z3,z4)
  u=mvrnorm(n,rep(0,4),covmatrix)
  w=z+u
  y=rep(0,n)
  error=rnorm(n,0,sigma)## sd 
  for(i in c(1:n)){
    y[i]=2*g[i,1]-g[i,2]+z[i,]%*%beta+error[i]+1.2*delta*(exp(sqrt(inpro(x[i,],x[i,],tm)))-1)
  }
  return(list(Y=y,X=x,Z=z,W=w))
}

#Case 1  1.2*delta*(exp(sqrt(inpro(x[i,],x[i,],tm)))-1)
#Case 2  delta*exp(abs(z4[i]-1))*2/3

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



  

#Benchmark method
tgridnum=100
tm=c(1:tgridnum)*0.01-0.005
alpha=2*sqrt(2)*sin(pi*tm/2)-sqrt(2)*sin(3*pi*tm/2)
beta=c(1,0.8,0.5,0.3)
B=300
cishu=500
sigma=0.5
m=2
testlevel=0.95
####testlevel=0.90
covmatrix=diag(c(0,0,0,0))

for (n in c(100)) {
  for (delta in c(0)) {
    timresult=c()
    for (tim in c(1:cishu)) {
      {
        data=generate(n,covmatrix,sigma,delta)
        Y=data$Y
        Z=data$Z
        X=data$X
        W=data$W
        H=pca(data$X,m,tm)$h
        eigf=pca(data$X,m,tm)$eig.functions
        S=H%*%solve(t(H)%*%H)%*%t(H)
        xnorm=c()
        for (i in c(1:n)) {
          xnorm[i]=sqrt(inpro(X[i,],X[i,],tm))
        }
        h1=sd(xnorm)*1.06*n^(-0.2)
        d=matrix(rep(0,n^2),n,n)
        for(i in c(2:n)){
          for(j in c(1:(i-1))){
            d[i,j]=sqrt(inpro(X[i,]-X[j,],X[i,]-X[j,],tm))
          }
        }
        hz24=sd(data$Z[,4])*1.06*n^(-0.2)
        hz23=sd(data$Z[,3])*1.06*n^(-0.2)
        hz22=sd(data$Z[,2])*1.06*n^(-0.2)
        hz21=sd(data$Z[,1])*1.06*n^(-0.2)
        Kz=matrix(rep(0,n^2),n,n)
        for(i in c(2:n)){
          for(j in c(1:(i-1))){
            # 4个分量都准确观测
            Kz[i,j]=fkernel(Z[i,4]-Z[j,4],hz24)*fkernel(d[i,j],h1)*fkernel(Z[i,1]-Z[j,1],hz21)*fkernel(Z[i,2]-Z[j,2],hz22)*fkernel(Z[i,3]-Z[j,3],hz23)
          }
        }
        Kz=Kz+t(Kz)
      }
      #BEestimation
      {
        betahat=solve(t(Z)%*%(diag(rep(1,n))-S)%*%Z)%*%t(Z)%*%(diag(rep(1,n))-S)%*%Y
        alphahat=solve(t(H)%*%H)%*%t(H)%*%(Y-Z%*%betahat)
        alphahatf=eigf%*%alphahat
        ehat=rep(0,n)
        for(i in c(1:n)){
          ehat[i]=Y[i]-t(Z[i,])%*%betahat-inpro(alphahatf,X[i,],tm)
        }
        statvalue=t(ehat)%*%Kz%*%ehat
      }
      ##BEbootstrap
      ehatcun=ehat
      bvalue=c()
      for(index in c(1:B)){
        eb=c()
        yb=c()
        for(i in c(1:n)){
          if(runif(1)<(5+sqrt(5))/10){
            eb[i]=(1-sqrt(5))/2
          }
          else
            eb[i]=(1+sqrt(5))/2
        }
        ebnew=eb*ehatcun
        for(i in c(1:n)){
          yb[i]=t(Z[i,])%*%betahat+inpro(alphahatf,X[i,],tm)+ebnew[i]
        }
        bbetahat=solve(t(Z)%*%(diag(rep(1,n))-S)%*%Z)%*%t(Z)%*%(diag(rep(1,n))-S)%*%yb
        balphahat=solve(t(H)%*%H)%*%t(H)%*%(yb-Z%*%bbetahat)
        balphahatf=eigf%*%balphahat
        ehat=rep(0,n)
        for(i in c(1:n)){
          ehat[i]=yb[i]-t(Z[i,])%*%bbetahat-inpro(balphahatf,X[i,],tm)
        }
        bvalue[index]=t(ehat)%*%Kz%*%ehat
      }
      timresult[tim]=(abs(statvalue)>quantile(abs(bvalue),testlevel))
      if(tim%%10==0){
        print(tim)
      }
    }
    print(c("Benckmark",covmatrix[1,1],n,delta,mean(timresult)))
  }
}
  
  
 
#Known method

for (n in c(100,200)) {
  for (ztvalue in c(0.8,0.4)) {
    covmatrix=diag(c(ztvalue,ztvalue,ztvalue,0))
    for (delta in c(0,1,1.5,2,2.5)) {
      timresult=c()
      for (tim in c(1:cishu)) {
        {
          data=generate(n,covmatrix,sigma,delta)
          Y=data$Y
          Z=data$Z
          X=data$X
          W=data$W
          H=pca(data$X,m,tm)$h
          eigf=pca(data$X,m,tm)$eig.functions
          S=H%*%solve(t(H)%*%H)%*%t(H)
          xnorm=c()
          for (i in c(1:n)) {
            xnorm[i]=sqrt(inpro(X[i,],X[i,],tm))
          }
          h1=sd(xnorm)*1.06*n^(-0.2)
          d=matrix(rep(0,n^2),n,n)
          for(i in c(2:n)){
            for(j in c(1:(i-1))){
              d[i,j]=sqrt(inpro(X[i,]-X[j,],X[i,]-X[j,],tm))
            }
          }
          hz24=sd(data$Z[,4])*1.06*n^(-0.2)
          Kz=matrix(rep(0,n^2),n,n)
          for(i in c(2:n)){
            for(j in c(1:(i-1))){
              # 只有第4个分量准确观测
              Kz[i,j]=fkernel(Z[i,4]-Z[j,4],hz24)*fkernel(d[i,j],h1)
            }
          }
          Kz=Kz+t(Kz)
        }
        #CPLSE1estimation
        {
          betahat=solve(t(W)%*%(diag(rep(1,n))-S)%*%W-n*covmatrix)%*%t(W)%*%(diag(rep(1,n))-S)%*%Y
          alphahat=solve(t(H)%*%H)%*%t(H)%*%(Y-W%*%betahat)
          alphahatf=eigf%*%alphahat
          ehat=rep(0,n)
          for(i in c(1:n)){
            ehat[i]=Y[i]-t(W[i,])%*%betahat-inpro(alphahatf,X[i,],tm)
          }
          statvalue=t(ehat)%*%Kz%*%ehat
        }
        ##CPLSE1bootstrap
        ehatcun=ehat
        bvalue=c()
        for(index in c(1:B)){
          eb=c()
          yb=c()
          for(i in c(1:n)){
            if(runif(1)<(5+sqrt(5))/10){
              eb[i]=(1-sqrt(5))/2
            }
            else
              eb[i]=(1+sqrt(5))/2
          }
          ebnew=eb*ehatcun
          for(i in c(1:n)){
            yb[i]=t(W[i,])%*%betahat+inpro(alphahatf,X[i,],tm)+ebnew[i]
          }
          bbetahat=solve(t(W)%*%(diag(rep(1,n))-S)%*%W)%*%t(W)%*%(diag(rep(1,n))-S)%*%yb
          balphahat=solve(t(H)%*%H)%*%t(H)%*%(yb-W%*%bbetahat)
          balphahatf=eigf%*%balphahat
          ehat=rep(0,n)
          for(i in c(1:n)){
            ehat[i]=yb[i]-t(W[i,])%*%bbetahat-inpro(balphahatf,X[i,],tm)
          }
          bvalue[index]=t(ehat)%*%Kz%*%ehat
        }
        timresult[tim]=(abs(statvalue)>quantile(abs(bvalue),testlevel))
      }
      print(c("Known",ztvalue,n,delta,mean(timresult)))
    }
  }
}
  



#naive method


for (n in c(100,200)) {
  for (ztvalue in c(0.8,0.4)) {
    covmatrix=diag(c(ztvalue,ztvalue,ztvalue,0))
    for (delta in c(0,1,1.5,2,2.5)) {
      timresult=c()
      for (tim in c(1:cishu)) {
        {
          data=generate(n,covmatrix,sigma,delta)
          Y=data$Y
          Z=data$Z
          X=data$X
          W=data$W
          H=pca(data$X,m,tm)$h
          eigf=pca(data$X,m,tm)$eig.functions
          S=H%*%solve(t(H)%*%H)%*%t(H)
          xnorm=c()
          for (i in c(1:n)) {
            xnorm[i]=sqrt(inpro(X[i,],X[i,],tm))
          }
          h1=sd(xnorm)*1.06*n^(-0.2)
          d=matrix(rep(0,n^2),n,n)
          for(i in c(2:n)){
            for(j in c(1:(i-1))){
              d[i,j]=sqrt(inpro(X[i,]-X[j,],X[i,]-X[j,],tm))
            }
          }
          hz24=sd(data$Z[,4])*1.06*n^(-0.2)
          Kz=matrix(rep(0,n^2),n,n)
          for(i in c(2:n)){
            for(j in c(1:(i-1))){
              Kz[i,j]=fkernel(Z[i,4]-Z[j,4],hz24)*fkernel(d[i,j],h1)
            }
          }
          Kz=Kz+t(Kz)
        }
        #naive estimation
        {
          betahat=solve(t(W)%*%(diag(rep(1,n))-S)%*%W)%*%t(W)%*%(diag(rep(1,n))-S)%*%Y
          alphahat=solve(t(H)%*%H)%*%t(H)%*%(Y-W%*%betahat)
          alphahatf=eigf%*%alphahat
          ehat=rep(0,n)
          for(i in c(1:n)){
            ehat[i]=Y[i]-t(W[i,])%*%betahat-inpro(alphahatf,X[i,],tm)
          }
          statvalue=t(ehat)%*%Kz%*%ehat
        }
        ##naive bootstrap
        ehatcun=ehat
        bvalue=c()
        for(index in c(1:B)){
          eb=c()
          yb=c()
          for(i in c(1:n)){
            if(runif(1)<(5+sqrt(5))/10){
              eb[i]=(1-sqrt(5))/2
            }
            else
              eb[i]=(1+sqrt(5))/2
          }
          ebnew=eb*ehatcun
          for(i in c(1:n)){
            yb[i]=t(W[i,])%*%betahat+inpro(alphahatf,X[i,],tm)+ebnew[i]
          }
          bbetahat=solve(t(W)%*%(diag(rep(1,n))-S)%*%W)%*%t(W)%*%(diag(rep(1,n))-S)%*%yb
          balphahat=solve(t(H)%*%H)%*%t(H)%*%(yb-W%*%bbetahat)
          balphahatf=eigf%*%balphahat
          ehat=rep(0,n)
          for(i in c(1:n)){
            ehat[i]=yb[i]-t(W[i,])%*%bbetahat-inpro(balphahatf,X[i,],tm)
          }
          bvalue[index]=t(ehat)%*%Kz%*%ehat
        }
        timresult[tim]=(abs(statvalue)>quantile(abs(bvalue),testlevel))
      }
      print(c("Naive",ztvalue,n,delta,mean(timresult)))
    }
  }
}
  
  
  
#UnKnown method

  
for (n in c(100,200)) {
  for (ztvalue in c(0.8,0.4)) {
    covmatrix=diag(c(ztvalue,ztvalue,ztvalue,0))
    for (delta in c(0,1,1.5,2,2.5)) {
      timresult=c()
      for (tim in c(1:cishu)) {
        {
          data=generate(n,covmatrix,sigma,delta)
          Y=data$Y
          Z=data$Z
          X=data$X
          W=data$W
          H=pca(data$X,m,tm)$h
          eigf=pca(data$X,m,tm)$eig.functions
          S=H%*%solve(t(H)%*%H)%*%t(H)
          xnorm=c()
          for (i in c(1:n)) {
            xnorm[i]=sqrt(inpro(X[i,],X[i,],tm))
          }
          h1=sd(xnorm)*1.06*n^(-0.2)
          d=matrix(rep(0,n^2),n,n)
          for(i in c(2:n)){
            for(j in c(1:(i-1))){
              d[i,j]=sqrt(inpro(X[i,]-X[j,],X[i,]-X[j,],tm))
            }
          }
          hz24=sd(data$Z[,4])*1.06*n^(-0.2)
          Kz=matrix(rep(0,n^2),n,n)
          for(i in c(2:n)){
            for(j in c(1:(i-1))){
              Kz[i,j]=fkernel(Z[i,4]-Z[j,4],hz24)*fkernel(d[i,j],h1)
            }
          }
          Kz=Kz+t(Kz)
        }
        #CPLSE2estimation
        {
          errva1=mvrnorm(n,rep(0,4),covmatrix)
          errva2=mvrnorm(n,rep(0,4),covmatrix)
          errva3=mvrnorm(n,rep(0,4),covmatrix)
          storage=matrix(rep(0,16),4,4)
          for(index in c(1:n)){
            storage=storage+(2*errva1[index,]-errva2[index,]-errva3[index,])%*%t(2*errva1[index,]-errva2[index,]-errva3[index,])
            storage=storage+(2*errva2[index,]-errva1[index,]-errva3[index,])%*%t(2*errva2[index,]-errva1[index,]-errva3[index,])
            storage=storage+(2*errva3[index,]-errva1[index,]-errva2[index,])%*%t(2*errva3[index,]-errva2[index,]-errva2[index,])
          }
          covmatrixhat=storage/n/18
          betahat=solve(t(W)%*%(diag(rep(1,n))-S)%*%W-n*covmatrixhat)%*%t(W)%*%(diag(rep(1,n))-S)%*%Y
          alphahat=solve(t(H)%*%H)%*%t(H)%*%(Y-W%*%betahat)
          alphahatf=eigf%*%alphahat
          ehat=rep(0,n)
          for(i in c(1:n)){
            ehat[i]=Y[i]-t(W[i,])%*%betahat-inpro(alphahatf,X[i,],tm)
          }
          statvalue=t(ehat)%*%Kz%*%ehat
        }
        ##CPLSE2bootstrap
        ehatcun=ehat
        bvalue=c()
        for(index in c(1:B)){
          eb=c()
          yb=c()
          for(i in c(1:n)){
            if(runif(1)<(5+sqrt(5))/10){
              eb[i]=(1-sqrt(5))/2
            }
            else
              eb[i]=(1+sqrt(5))/2
          }
          ebnew=eb*ehatcun
          for(i in c(1:n)){
            yb[i]=t(W[i,])%*%betahat+inpro(alphahatf,X[i,],tm)+ebnew[i]
          }
          bbetahat=solve(t(W)%*%(diag(rep(1,n))-S)%*%W)%*%t(W)%*%(diag(rep(1,n))-S)%*%yb
          balphahat=solve(t(H)%*%H)%*%t(H)%*%(yb-W%*%bbetahat)
          balphahatf=eigf%*%balphahat
          ehat=rep(0,n)
          for(i in c(1:n)){
            ehat[i]=yb[i]-t(W[i,])%*%bbetahat-inpro(balphahatf,X[i,],tm)
          }
          bvalue[index]=t(ehat)%*%Kz%*%ehat
        }
        timresult[tim]=(abs(statvalue)>quantile(abs(bvalue),testlevel))
      }
      print(c("UnKnown",covmatrix[1,1],n,delta,mean(timresult)))
    } 
  }
}
  
  