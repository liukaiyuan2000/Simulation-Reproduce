rm(list=ls())
set.seed(123456)
library(MASS)

 S1 <- function(data1,theta){
   # theta <- c(beta0,1)
   #theta <- theta.old
   X=data1$X
   Y=data1$Y
   delta=data1$delta
   p <- ncol(X)
   n <- nrow(X)
   eps <- 1e-3/n
   index <- (Y-X%*%theta[1:p])/(theta[p+1]+1e-3/n)
   temp <- delta*(index/theta[p+1])+(1-delta)*dnorm(index)/(pnorm(-index)+eps)
   SS1 <- apply(t(X)%*%diag(as.vector(temp)),1,sum) 
   SS2 <- sum(delta*(-1/theta[p+1]+index^2/theta[p+1])+(1-delta)*dnorm(index)*index/theta[p+1]/(pnorm(-index)+eps))
   c(SS1,SS2)
   
 }
 
NH <- function(data1,theta){
   X=data1$X
   Y=data1$Y
   delta=data1$delta
   p <- ncol(X)
   n <- nrow(X)
   eps <- 1e-3/n
   index <- (Y-X%*%theta[1:p])/(theta[p+1]+1e-3/n)
   temp1 <- -t(X)%*%diag(as.vector(delta/theta[p+1]/theta[p+1]))%*%X
   Temp <- -(1-delta)*index*dnorm(index)/(pnorm(-index)+eps)/theta[p+1]+pnorm(index)^2/theta[p+1]
   temp2 <- -t(X)%*%diag(as.vector(Temp))%*%X
   H11 <- temp1+temp2
   temp1 <- delta*(1-3*index^2)/theta[p+1]/theta[p+1]
   temp21 <- (1-delta)*(index*dnorm(index)*index^2-2*dnorm(index)*index)/theta[p+1]/theta[p+1]
   temp22 <- -(1-delta)*dnorm(index)^2*index^2/theta[p+1]/theta[p+1]/((pnorm(-index)+eps)^2)
   
   H22 <- sum(temp1)+sum(temp21/(pnorm(-index)+eps))+sum(temp22)
   temp1 <- -t(X)%*%(index*delta/theta[p+1])
   temp21 <- t(X)%*%(index*dnorm(index)*index/(pnorm(-index)+eps))
   temp22 <- -t(X)%*%((dnorm(index))^2*index^2/theta[p+1]/theta[p+1]/((pnorm(-index)+eps)^2))
   
   H21 <- temp21+temp22
   
   rbind(cbind(H11,(H21)),c(H21,H22))
   
}
 
n <- 200
beta0 <- c(1,-1)
 
p <- length(beta0)
Sigma <- matrix(0,p,p)
Sigma <- 0.5^(abs(col(Sigma)-row(Sigma)))
Sim <- 10
THETA <- vector("list",4)
for(sim in 1:Sim){
   cat(sim,"\r")
   X <- mvrnorm(n,rep(0,p),Sigma)
   Ystar <- X%*%beta0+rnorm(n)
   C <- runif(n,0,2)
   Y <- pmin(C,Ystar)
   delta <- 1*(Ystar<C)
   mean(delta)
   data1 =list(X=X,Y=Y,delta=delta)
   
   fit1 <- lm(Y~X+0)
   # summary(fit1)
   CY <- Y[delta==1]
   CX <- X[delta==1,]
   fit2 <- lm(CY~CX+0)
   beta.naive2 <- fit2$coef
   beta.naive1 <- fit1$coef
   sigma.naive1 <- sqrt(sum(fit1$residuals^2)/(n-p))
   sigma.naive2 <- sqrt(sum(fit2$residuals^2)/(n-p))
   theta.naive1 <- c(beta.naive1,sigma.naive1^2)
   theta.naive2 <- c(beta.naive2,sigma.naive2^2)
   
   fit <- lm(Ystar~X+0)
   
   beta.gold <- fit$coef
   
   
   sigma.gold <- sqrt(sum(fit$residuals^2)/(n-p))
   theta.gold <- c(beta.gold,sigma.gold^2)
   
   theta.new <- c(lm(CY~CX+0)$coef,1)
   
   theta.old <-  runif(p+1)
   Loop <- 10
   loop <- 0
   er <- 1
   while(loop<Loop&er>1e-5){
     loop <- loop+1
     
     theta.old <- theta.new
     #theta.old <- c(beta0,1) 
     SS <- S1(data1,theta.old)
     H <- NH(data1,theta.old)
     #theta <- c(beta0,1)
     #NH(data1,c(beta0,1))/n
     theta.new <- theta.old-solve(H)%*%SS#更新方程
     
     er <- sum(abs(theta.old-theta.new))
     #cat(c(theta.new,loop,er),"\n")
   }
   
   #theta.new <- as.vector(theta.new)
   colnames(theta.new) <- "theta.new"
   THETA[[1]] <- rbind(THETA[[1]],theta.naive1-c(beta0,1))
   THETA[[2]] <- rbind(THETA[[2]],theta.naive2-c(beta0,1))
   THETA[[3]] <- rbind(THETA[[3]],t(theta.new)-c(beta0,1))
   THETA[[4]] <- rbind(THETA[[4]],theta.gold-c(beta0,1))
   
   
 }

c1 <- c(apply(THETA[[1]],2,mean),apply(THETA[[1]],2,sd))
c2 <- c(apply(THETA[[2]],2,mean),apply(THETA[[2]],2,sd))
c3 <- c(apply(THETA[[3]],2,mean),apply(THETA[[3]],2,sd))
c4 <- c(apply(THETA[[4]],2,mean),apply(THETA[[4]],2,sd))





res <- rbind( c1,c2,c3,c4)
print(res)

 
 
 
 