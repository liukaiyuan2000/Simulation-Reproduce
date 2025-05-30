
############################
#computational statistics 
##Exercise 2.3##
library(MASS)
n <- 500
alpha0 <- 2
beta0 <- c(-1,1)
p <- length(beta0)
Sigmma <- matrix(0,p,p)
Sigmma <- 0.5^abs(row(Sigmma)-col(Sigmma)) 
X <- mvrnorm(n, rep(0,p),Sigmma)
X[,1] <- 1*(X[,1]<0)
index <- X%*%beta0
U <- runif(n)
Ystar <- -log(U)*exp(-index)
Ystar <- Ystar^(1/alpha0)
C <- runif(n,0,7)
delta <- 1*(Ystar<C)
Y <- pmin(Ystar,C) 

data1=list(Y=Y,delta=delta,X=X)
S1 <- function(data1,theta){
  Y=as.vector(data1$Y)
  X=data1$X
  delta=data1$delta
  index <-  X%*%theta[-1]
  c1=sum(-Y^theta[1]*log(Y)*exp(-index)+delta/theta[1]+delta*log(Y))
  c2=apply(diag(as.vector(Y^theta[1]*exp(-index)-delta))%*%X,2,sum)
  c(c1,c2) 
}
H1 <- function(data1,theta){
  Y=data1$Y
  X=data1$X
  delta=data1$delta
  index <- X%*%theta[-1]
  H12 <- apply(diag(as.vector(Y^theta[1]*log(Y)*exp(-index)))%*%X,2,sum)
  H11 <- sum(-Y^theta[1]*(log(Y))^2*exp(-index)-delta/theta[1]/theta[1]) 
  H22 <- t(X)%*%diag(as.vector(-Y^theta[1]*exp(-index)))%*%X
  rbind(c(H11,H12),cbind((H12),H22))
}


#theta <- theta.old

eps <- 10
loop <- 0
Loop <- 1000
theta.new <- runif(3)
theta.old <- theta.new
# theta <- c(2,  1,-1)
while(loop<Loop&eps>1e-3){
  loop <- loop+1
  theta.old <- theta.new
  S <-S1(data1,theta.old) 
  H <-H1(data1,theta.old) 
  theta.new <- as.vector(theta.old)-solve(H)%*%S
  theta.new <- as.vector( theta.new)
  eps <- sum((theta.old - theta.new)^2)
  eps <- sqrt(eps)
  cat(c(eps,theta.new),"\n")
}
print(theta.new)

####################################
rm(list=ls())
#Y=min(T,C),delta=1*(T<=C)
Yt=c(6,6,6,6,7,9,10,10,11,13,16,17,19,20,22,23,25,32,32,34,35)
Yc=c(1,1,2,2,3,4,4,5,5,8,8,8,8,11,11,12,12,15,17,22,23)
delta <- c(0,1,1,1,1,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,rep(1,length(Yc)))
Y <- c(Yt,Yc)
X <- c(rep(1,length(Yt)),rep(0,length(Yc)))
X <- cbind(1,X)
data1=list(Y=Y,delta=delta,X=X)
S1 <- function(data1,theta){
  Y=as.vector(data1$Y)
  X=data1$X
  delta=data1$delta
  index <-  X%*%theta[-1]
  c1=sum(-Y^theta[1]*log(Y)*exp(-index)+delta/theta[1]+delta*log(Y))
  c2=apply(diag(as.vector(Y^theta[1]*exp(-index)-delta))%*%X,2,sum)
  c(c1,c2) 
}
H1 <- function(data1,theta){
  Y=data1$Y
  X=data1$X
  delta=data1$delta
  index <- X%*%theta[-1]
  H12 <- apply(diag(as.vector(Y^theta[1]*log(Y)*exp(-index)))%*%X,2,sum)
  H11 <- sum(-Y^theta[1]*(log(Y))^2*exp(-index)-delta/theta[1]/theta[1]) 
  H22 <- t(X)%*%diag(as.vector(-Y^theta[1]*exp(-index)))%*%X
  rbind(c(H11,H12),cbind((H12),H22))
}


#theta <- theta.old

eps <- 10
loop <- 0
Loop <- 100

theta.new <-c(1,1,1)
theta.old <- theta.new

while(loop<Loop&eps>1e-5){
  loop <- loop+1
  theta.old <- theta.new
  
  S <-S1(data1,theta.old)
  H <-H1(data1,theta.old)
  theta.new <- as.vector(theta.old)-solve(H)%*%S
  theta.new <- as.vector( theta.new)
  eps <- sum((theta.old - theta.new)^2)
  eps <- sqrt(eps)
  #cat(c(eps,theta.new),"\n")
}
print(round(theta.new,3))

cov1=solve(-H)
theta.sd <- sqrt(diag(cov1))
Tn <- theta.new/theta.sd
pvalue <-2*(1-pnorm(Tn))
cbind(theta.new,theta.sd,pvalue)



p <- length(theta.new)
rcov <- matrix(0,p,p)
for(i in 1:p)
  for(j in 1:p)
    rcov[i,j]=cov1[i,j]/sqrt(cov1[i,i])/sqrt(cov1[j,j])





