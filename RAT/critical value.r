# A rank-based adaptive independence test for high-dimensional data
rm(list=ls())
setwd("D:/Desktop/RAT")
Rcpp::sourceCpp("allcpp.cpp")
set.seed(123456)  
library(RcppArmadillo)
library(Rcpp)
library(MASS)  
t1 <- Sys.time()
RES <- NULL
N <- 200
P <- 400
 
Sim <- 2000
alpha <- 0.05
ECPRES <- NULL
cv <- NULL
for(n in N){
  REs <- NULL
  for(p in P){
    ECP <- matrix(0,Sim,7)
    Pnp<-rep(0,6)
    Tn <- rep(0,6)
    for(sim in 1:Sim){
      cat(n,p,sim,"\r")
      X <- matrix(runif(n*p),n,p) 
      Tn <- allmain(X)
      #sum
      ECP[sim,1:3]<-Tn[1:3]
      #max
      ECP[sim,4:6]<-Tn[4:6]-4*log(p)+log(log(p))
      #sum-sum
      Pnp[1] <- 2*(1-pnorm(abs(Tn[1])))
      Pnp[2] <- 2*(1-pnorm(abs(Tn[2])))
      Pnp[3] <- 2*(1-pnorm(abs(Tn[3])))
      Pnp[4] <- 1-exp(-exp(-(Tn[4]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
      Pnp[5] <- 1-exp(-exp(-(Tn[5]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
      Pnp[6] <- 1-exp(-exp(-(Tn[6]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
      Cnp=min(Pnp[1],Pnp[2],Pnp[3],Pnp[4],Pnp[5],Pnp[6])
      ECP[sim,7]<-Cnp
  
    }
    temp1 <- c(quantile(ECP[,1], c(alpha/2,1-alpha/2)),quantile(ECP[,2], c(alpha/2,1-alpha/2)),quantile(ECP[,3], c(alpha/2,1-alpha/2)))
    temp2 <- c(quantile(ECP[,4],1-alpha),quantile(ECP[,5],1-alpha),quantile(ECP[,6],1-alpha),quantile(ECP[,7],alpha))
    
    
    
    cv <- rbind(cv,c(temp1,temp2))
  }
  
}
#write(round(ECPRES,3),"Tnstartest.txt")
#write(cv,"Tnstarcvn.txt")
t2 <- Sys.time()
print((t2-t1))
write.table(cv, file = "cv2.txt")

#write.csv(cv,file = "C:/Users/Administrator/Documents/Xiangyu Shi/CSSC/cv1.csv",row.names = FALSE)



