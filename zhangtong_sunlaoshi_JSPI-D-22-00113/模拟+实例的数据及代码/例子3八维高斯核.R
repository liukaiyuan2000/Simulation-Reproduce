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
    y[i]=g[i,1]+g[i,2]+z[i,]%*%beta+error[i]+delta*(z[i,8]^3)/4
  }
  return(list(Y=y,X=x,Z=z,W=w))
}

#different choices of the deviation function
#Case 1 delta*4*log(inpro(x[i,],x[i,],tm)+1)
#Case 2 delta*(z[i,8]^3)/4

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



{
#Benchmark method
tgridnum=100
tm=c(1:tgridnum)*0.01-0.005
#alpha=sqrt(2)*sin(pi*tm/2)+sqrt(2)*sin(3*pi*tm/2)这行没用到
beta=c(0.8,0.6,0.4,0.2,0.8,0.6,0.4,0.2)
B=300
cishu=500
sigma=0.5
m=2
covmatrix=diag(c(0.8,0.8,0.8,0.8,0.8,0.8,0,0))
testlevel=0.90
####testlevel=0.90

for (n in c(200)) {
  for (delta in c(1,1.5)) {
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
        # hz1=sd(data$Z[,1])*1.06*n^(-0.2)
        # hz2=sd(data$Z[,2])*1.06*n^(-0.2)
        # hz3=sd(data$Z[,3])*1.06*n^(-0.2)
        # hz4=sd(data$Z[,4])*1.06*n^(-0.2)
        # hz5=sd(data$Z[,5])*1.06*n^(-0.2)
        # hz6=sd(data$Z[,6])*1.06*n^(-0.2)
        hz7=sd(data$Z[,7])*1.06*n^(-0.2)
        hz8=sd(data$Z[,8])*1.06*n^(-0.2)
        Kz=matrix(rep(0,n^2),n,n)
        for(i in c(2:n)){
          for(j in c(1:(i-1))){
            Kz[i,j]=fkernel(d[i,j],h1)*fkernel(Z[i,7]-Z[j,7],hz7)*fkernel(Z[i,8]-Z[j,8],hz8)#*fkernel(Z[i,1]-Z[j,1],hz1)*fkernel(Z[i,2]-Z[j,2],hz2)*fkernel(Z[i,3]-Z[j,3],hz3)*fkernel(Z[i,4]-Z[j,4],hz4)*fkernel(Z[i,5]-Z[j,5],hz5)*fkernel(Z[i,6]-Z[j,6],hz6)
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
      # if(tim%%50==0){
      #   print(tim)
      # }
    }
    print(c("Benckmark",covmatrix[1,1],n,delta,mean(timresult)))
  }
}



#Known method
####ztvalue means the covariance of the measurement error
for (n in c(200)){
  for (ztvalue in c(0.4)) {
    covmatrix=diag(c(ztvalue,ztvalue,ztvalue,ztvalue,ztvalue,ztvalue,0,0))
    for (delta in c(0,0,0,1.5,2,2.5)) {
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
          hz7=sd(data$Z[,7])*1.06*n^(-0.2)
          hz8=sd(data$Z[,8])*1.06*n^(-0.2)
          Kz=matrix(rep(0,n^2),n,n)
          for(i in c(2:n)){
            for(j in c(1:(i-1))){
              Kz[i,j]=fkernel(d[i,j],h1)*fkernel(Z[i,7]-Z[j,7],hz7)*fkernel(Z[i,8]-Z[j,8],hz8)
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
        # if(tim%%100==0){
        #   print(tim)
        # }
      }
      print(c("Known",covmatrix[1,1],n,delta,mean(timresult)))
    }
  }
}


#naive method
for (n in c(100,200)) {
  for (ztvalue in c(0.8,0.4)) {
    covmatrix=diag(c(ztvalue,ztvalue,ztvalue,ztvalue,ztvalue,ztvalue,0,0))
    for (delta in c(1,1.5,2,2.5)) {
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
          ##########DA XIU
          # hz1=sd(data$W[,1])*1.06*n^(-0.2)
          # hz2=sd(data$W[,2])*1.06*n^(-0.2)
          # hz3=sd(data$W[,3])*1.06*n^(-0.2)
          # hz4=sd(data$W[,4])*1.06*n^(-0.2)
          # hz5=sd(data$W[,5])*1.06*n^(-0.2)
          # hz6=sd(data$W[,6])*1.06*n^(-0.2)
          hz7=sd(data$W[,7])*1.06*n^(-0.2)
          hz8=sd(data$W[,8])*1.06*n^(-0.2)
          Kz=matrix(rep(0,n^2),n,n)
          for(i in c(2:n)){
            for(j in c(1:(i-1))){
              Kz[i,j]=fkernel(d[i,j],h1)*fkernel(W[i,7]-W[j,7],hz7)*fkernel(W[i,8]-W[j,8],hz8)#*fkernel(W[i,1]-W[j,1],hz1)*fkernel(W[i,2]-W[j,2],hz2)*fkernel(W[i,3]-W[j,3],hz3)*fkernel(W[i,4]-W[j,4],hz4)*fkernel(W[i,5]-W[j,5],hz5)*fkernel(W[i,6]-W[j,6],hz6)
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
         # if(tim%%50==0){
         #   print(tim)
         # }
      }
      print(c("Naive",covmatrix[1,1],n,delta,mean(timresult)))
    }
  }
}




#UnKnown method
for (n in c(200)) {
  for (ztvalue in c(0.4)) {
    covmatrix=diag(c(ztvalue,ztvalue,ztvalue,ztvalue,ztvalue,ztvalue,0,0))
    for (delta in c(0,0,0,0,0,0,0)) {
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
          hz7=sd(data$Z[,7])*1.06*n^(-0.2)
          hz8=sd(data$Z[,8])*1.06*n^(-0.2)
          Kz=matrix(rep(0,n^2),n,n)
          for(i in c(2:n)){
            for(j in c(1:(i-1))){
              Kz[i,j]=fkernel(d[i,j],h1)*fkernel(Z[i,7]-Z[j,7],hz7)*fkernel(Z[i,8]-Z[j,8],hz8)
            }
          }
          Kz=Kz+t(Kz)
        }
        #CPLSE1estimation
        {
          errva1=mvrnorm(n,rep(0,8),covmatrix)
          errva2=mvrnorm(n,rep(0,8),covmatrix)
          errva3=mvrnorm(n,rep(0,8),covmatrix)
          storage=matrix(rep(0,64),8,8)
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
        if(tim%%50==0){
          print(tim)
        }
      }
      print(c("UnKnown",covmatrix[1,1],n,delta,mean(timresult)))
    } 
  }
}





