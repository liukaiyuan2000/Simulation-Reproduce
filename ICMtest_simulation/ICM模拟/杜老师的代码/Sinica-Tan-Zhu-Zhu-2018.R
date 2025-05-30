rm(list=ls())
library(MASS)
setwd('D:/桌面/ICM模拟/杜老师的代码')
Rcpp::sourceCpp("GWZtest.cpp")
Rcpp::sourceCpp("TZZtest.cpp")
Rcpp::sourceCpp("SGPtest.cpp")
t1 <- Sys.time()
N <- c(50,100)
A <- seq(0,0.4,by=0.2)
Sim <- 1000 
p <- 8
beta0 <- rep(1,p)/sqrt(p)
Sigma <- matrix(0,p,p)
Sigma <- 0.5^abs(col(Sigma)-row(Sigma))
RES <- NULL
RHO <- c(0,0.5)
for(rho in RHO){
  Sigma <- rho^abs(col(Sigma)-row(Sigma))
  REs <- NULL
  for(n in N){
    res=NULL
    for(a in A){
      #cat(rho,n,a,"\r")
      ECP <- matrix(0,Sim,3)
      for(sim in 1:Sim){
        cat(rho,n,a,sim,"\r")
        #X <- matrix(rnorm(n*p),n,p) # sigma=I
        X <- mvrnorm(n, rep(0,p), Sigma) # sigma = Sigma
        index <- X%*%beta0
        e <- rnorm(n)
        #Y <- index+a*cos(index*pi/2)+e #H21
        Y <- index+a*exp(index)/4+e #H22
        fit0 <- GWZtest(Y,X)
        fit1 <- TZZtest(Y,X) 
        fit2 <- SGPtest(Y,X)
        ECP[sim,] <- c(1*(abs(fit0)>qnorm(0.975)),1*(fit1[2]<0.05),1*(fit2[2]<0.05))
      }
      res=rbind(res,colMeans(ECP))
    }
    REs <- cbind(REs,res)
  }
  RES <- rbind(RES,REs)
  
}
 
RES=RES[,c(1,4,2,5,3,6)]
colnames(RES) <- c("GWZn=50","GWZn=100","TZZn=50","TZZn=100","SGPn=50","SGPn=100")
t2 <- Sys.time()
print(round(RES,4))
print(t2-t1)
