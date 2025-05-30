
Genedata <- function(n,beta0,rho=0.5){
  beta0 <- as.vector(beta0)
  p <- length(beta0)
  Sigmma <- matrix(0,p,p)
  Sigmma <- (rho)^abs(row(Sigmma)-col(Sigmma)) 
  X <- mvrnorm(n, rep(0,p),Sigmma)
 
  Y=X%*%beta0+rnorm(n)
 
  out=list(X=X,Y=Y)
  out
}
softfunction <- function(z,u){
  max(0,abs(z)-u)*sign(z)
  
}

CDLASSO <- function(mydata,lambda=NULL,eps=1e-5,Loop=1000){
  #Variable selection for linear regression model via the cyclic coordinate descent algorithm 
  #nonlinear jacobi iteration
  X <- mydata$X
  Y <- mydata$Y
  fit <- lm(Y~X+0)
  beta.new <- fit$coef
  beta.ls <- beta.new
  beta.old <- beta.new
  if(is.null(lambda)) lambda <-   var(fit$res)*log(n)/2
  er <- 1
  loop <- 0
  p <- ncol(X)
  while(er>eps&loop<Loop){
    loop <- loop+1
    beta.old <- beta.new
    r <- Y-X%*%beta.old
    for(i in 1:p){
      #r <- Y-X%*%beta.new
      temp <- drop(sum(X[,i]*(r+X[,i]*beta.old[i])))
      beta.new[i] <- softfunction(temp,lambda/(abs(beta.old[i])+1/n/n) ) /sum(X[,i]^2)
    }
    er <- sqrt(sum((beta.new-beta.old)^2))
  }
  df=sum(beta.new!=0)
  sigmma2.hat <- sum((Y-X%*%beta.new)^2)/(n-df)
 out=list(beta.hat=beta.new,sigmma2.hat=sigmma2.hat,beta.ls=beta.ls,loop=loop,er=er)
 out
  
}

CDLASSO1 <- function(mydata,lambda=NULL,eps=1e-5,Loop=1000){
  #Variable selection for linear regression model via the cyclic coordinate descent algorithm 
  #nonlinear Gauss-Seidel iteration
  X <- mydata$X
  Y <- mydata$Y
  fit <- lm(Y~X+0)
  beta.new <- fit$coef
  beta.ls <- beta.new
  beta.old <- beta.new
  if(is.null(lambda)) lambda <-   var(fit$res)*log(n)/2
  er <- 1
  loop <- 0
  p <- ncol(X)
  while(er>eps&loop<Loop){
    loop <- loop+1
    beta.old <- beta.new
    r <- Y-X%*%beta.old
    for(i in 1:p){
      #r <- Y-X%*%beta.new
      temp <- drop(sum(X[,i]*(r+X[,i]*beta.old[i])))
      beta.new[i] <- softfunction(temp,lambda/(abs(beta.old[i])+1/n/n) ) /sum(X[,i]^2)
    }
    er <- sqrt(sum((beta.new-beta.old)^2))
  }
  df=sum(beta.new!=0)
  sigmma2.hat <- sum((Y-X%*%beta.new)^2)/(n-df)
  out=list(beta.hat=beta.new,sigmma2.hat=sigmma2.hat,beta.ls=beta.ls,loop=loop,er=er)
  out
  
}


library(MASS)
n <- 200
beta0 <- c(1,0,0,0,0,0,2.2,0,0,0.-3.5,0,0,0,1.75)
mydata <- Genedata(n,beta0)

fit1 <- CDLASSO(mydata)
fit2 <- CDLASSO1(mydata)


rbind(fit1$beta.hat,fit2$beta.hat,fit1$beta.ls)

 

