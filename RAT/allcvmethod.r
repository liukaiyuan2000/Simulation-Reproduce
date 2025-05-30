# A rank-based adaptive independence test for high-dimensional data
rm(list=ls())
set.seed(123456)
library(MASS)
setwd("C:/Users/Administrator/Documents/Xiangyu Shi/CSSC")  
Rcpp::sourceCpp("allcpp.cpp")
t1 <- Sys.time()
RES <- NULL
N <- c(50,100,200)
P <- c(50,100,200,400)
Sim <- 2000 
alpha <- 0.05
cv <- read.table("cv2.txt")
round(cv,5)
cv <- as.matrix(cv)
m=0
for(n in N){
  REs <- NULL
  for(p in P){ 
    
    ECP <- matrix(0,Sim,7)
    Pnp<-rep(0,6)
    Tn <- rep(0,6)
    m=m+1
    cvs <- cv[m,]
    cvs <- drop(cvs)
    for(sim in 1:Sim){
      cat(n,p,sim,"\r")
      #H0
      #??????̬?ֲ???????
      X <- matrix(rnorm(n*p),n,p) 
      
      
      #?????ֲ???????
       #X <- matrix(rcauchy(n*p,0,1),n,p) 
      
      
      #???Ϸֲ???????
      # X1 <- rnorm(n*p)
      # index <- sample(1:(n*p),0.05*n*p)
      # X2<- rnorm(0.05*n*p,0,10)
      # X1[index] <- X1[index]+X2
      # X <- matrix(X1,n,p)
      
      
      #dense H1
      # ??????̬?ֲ???????
      #rho=0.03
      #Sigmma <- (1-rho)*diag(p)+rho*matrix(1,p,p)
      #X <- mvrnorm(n, rep(0,p),Sigmma)
      
      #???ɿ????ֲ???????
      # Y<-matrix(0,p,n)
      # Z<-matrix(rcauchy(n*p,0,1),n,p)
      # for (k in 1:p) {
      #   Y[k,]<-rowSums(Z[,-k])
      # }
      # X<-Z+t(Y)/10/p
      
      
      
      # ?????쳣�???????
      # rho=0.03
      # Sigmma <- (1-rho)*diag(p)+rho*matrix(1,p,p)
      # X1 <- mvrnorm(n, rep(0,p),Sigmma)
      # X1<-c(X1)
      # index <- sample(1:(n*p),0.05*n*p)
      # X2<- rnorm(0.05*n*p,0,10)
      # X1[index] <- X1[index]+X2
      # X <- matrix(X1,n,p)
      
      
      #sparse H1
      #??????̬?ֲ???????
      # sigmma=matrix(0,p,p)
      # 
      # Sigmma <- 1*diag(p)
      # Sigmma[2,1]<-2.5*sqrt(log(p)/n)
      # Sigmma[1,2]<-2.5*sqrt(log(p)/n)
      # 
      # X <- mvrnorm(n, rep(0,p),Sigmma)
      
      #???ɿ????ֲ???????
      # Z=matrix(rcauchy(n*p,0,1),n,p)
      # X=cbind(Z[,1]+sqrt(log(p)/n)*Z[,2],Z[,2]+sqrt(log(p)/n)*Z[,1],Z[,(3:p)])
      
      
      
      #?????쳣�???????
      # sigmma=matrix(0,p,p)
      # 
      # Sigmma <- 1*diag(p)
      # Sigmma[2,1]<-2.5*sqrt(log(p)/n)
      # Sigmma[1,2]<-2.5*sqrt(log(p)/n)
      # X1 <- mvrnorm(n, rep(0,p),Sigmma)
      # X1<-c(X1)
      # index <- sample(1:(n*p),0.05*n*p)
      # X2<- rnorm(0.05*n*p,0,10)
      # X1[index] <- X1[index]+X2
      # X <- matrix(X1,n,p)
      
      Tn <- allmain(X)
      Tn <- drop(Tn)
      #sum
      ECP[sim,1]<-1*((Tn[1]<cvs[1])||( Tn[1]>cvs[2]))
      ECP[sim,2]<-1*((Tn[2]<cvs[3])||( Tn[2]>cvs[4]))
      ECP[sim,3]<-1*((Tn[3]<cvs[5])||( Tn[3]>cvs[6]))
      #max
      ECP[sim,4]<-1*((Tn[4]-4*log(p)+log(log(p)))>cvs[7])
      ECP[sim,5]<-1*((Tn[5]-4*log(p)+log(log(p)))>cvs[8])
      ECP[sim,6]<-1*((Tn[6]-4*log(p)+log(log(p)))>cvs[9])
      #sum-sum
      Pnp[1] <- 2*(1-pnorm(abs(Tn[1])))
      Pnp[2] <- 2*(1-pnorm(abs(Tn[2])))
      Pnp[3] <- 2*(1-pnorm(abs(Tn[3])))
      Pnp[4] <- 1-exp(-exp(-(Tn[4]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
      Pnp[5] <- 1-exp(-exp(-(Tn[5]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
      Pnp[6] <- 1-exp(-exp(-(Tn[6]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
      Cnp=min(Pnp[1],Pnp[2],Pnp[3],Pnp[4],Pnp[5],Pnp[6])
      ECP[sim,7]<-1*(Cnp<cvs[10])
      #all
      ECP[sim,]<-c(ECP[sim,1],ECP[sim,2],ECP[sim,3],ECP[sim,4],ECP[sim,5],ECP[sim,6],ECP[sim,7])
    }
    REs <- rbind(REs,colMeans(ECP))
  }
  rownames(REs) <- paste("p=",P)
  RES <- rbind(RES,REs)
}
colnames(RES) <- c("sumrho","sumtau","sumfoot","maxrho","maxtau","maxfoot","summax")

t2 <- Sys.time()
print((t2-t1))
print(RES)

write.csv(RES,file = "C:/Users/Administrator/Documents/Xiangyu Shi/CSSC/dnorm1.csv")



 
