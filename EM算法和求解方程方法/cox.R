library(survival)
data(rats)
data2 <- rats
data2$sex[data2$sex=="f"] <- 1
data2$sex[data2$sex=="m"] <- 2
mod <- coxph(Surv(time,status)~.,data = data2)
summary(mod)
Y <- data2$time
X <- cbind(as.numeric(data2$rx),as.numeric(data2$sex))
delta <- data2$status
data1=list(Y=Y,delta=delta,X=X)
S1 <- function(data1,theta){
  Y=as.vector(data1$Y)
  X=data1$X
  delta=data1$delta
  eta <- 0
  for(i in 1:length(Y)){
    la <- 0
    for(l in 1:length(Y)){
      la <- la+(Y[i]<=Y[l])*(delta[l]==1)*as.vector(exp(t(X[l,])%*%theta))
    }
    eta[i] <- la
  }
  c <- rep(0,ncol(X))
  for(i in 1:length(Y)){
    c1 <- rep(0,ncol(X))
    for(l in 1:length(Y)){
      c1 <- c1+(Y[i]<=Y[l])*(delta[l]==1)*as.vector(exp(t(X[l,])%*%theta))*X[l,]
    }
    c <- c-delta[i]*(c1/eta[i]+X[i,])
  }
  return(c)
}
H1 <- function(data1,theta){
  Y=as.vector(data1$Y)
  X=data1$X
  delta=data1$delta
  eta <- 0
  for(i in 1:length(Y)){
    la <- 0
    for(l in 1:length(Y)){
      la <- la+(Y[i]<=Y[l])*(delta[l]==1)*as.vector(exp(X[l,]%*%theta))
    }
    eta[i] <- la
  }
  eta1 <- matrix(0,nrow=length(Y),ncol=ncol(X))
  for(i in 1:length(Y)){
    la <- rep(0,ncol(X))
    for(l in 1:length(Y)){
      la <- la+(Y[i]<=Y[l])*(delta[l]==1)*as.vector(exp(X[l,]%*%theta))*X[l,]
    }
    eta1[i,] <- la
  }
  H <- matrix(0,ncol(X),ncol(X))
  for(i in 1:length(Y)){
    H1 <- matrix(0,ncol(X),ncol(X))
    for(l in 1:length(Y)){
      H1 <- H1+(Y[i]<=Y[l])*(delta[l]==1)*as.vector(exp(X[l,]%*%theta))*X[l,]%*%t(X[l,])
    }
    H <- H-delta[i]*(H1*eta[i]-eta1[i,]%*%t(eta1[i,]))/eta[i]^2
  }
  return(H)
}




eps <- 10
loop <- 0
Loop <- 50
theta.new <- c(1,1)
theta.old <- theta.new
while(loop<Loop&eps>1e-5){
  loop <- loop+1
  theta.old <- theta.new
  S <- S1(data1,theta.old)
  H <- H1(data1,theta.old)
  theta.new <- as.vector(theta.old)-solve(H)%*%S
  theta.new <- as.vector( theta.new)
  eps <- sum((theta.old - theta.new)^2)
  eps <- sqrt(eps)
}
#out1=list(theta.new=theta.new,loop=loop)
round(theta.new,5)

aa <- seq(-1000,0,by=1)
bb <- 0
for(i in 1:length(aa)){
  bb[i] <- S1(data1,aa[i])
}
plot(aa,bb)

eta <- 0
for(i in 1:length(Y)){
  la <- 0
  for(l in 1:length(Y)){
    la <- la+(Y[i]<=Y[l])*(delta[l]==1)*as.vector(exp(t(X[l,])%*%theta))
  }
  eta[i] <- la
}