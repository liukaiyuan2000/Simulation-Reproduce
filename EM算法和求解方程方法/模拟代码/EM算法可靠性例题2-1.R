mydata <- read.csv("¿É¿¿ÐÔ2.1.csv", header = T)
r <- length(delta)-length(delta[delta == 0])
n <- nrow(mydata)
st <- mydata$t
Y <- log(mydata$t)
X <- as.matrix(mydata[,2:8])
X[,1] <- X[,1]-mean(X[,1])
X[,2] <- X[,2]-mean(X[,2])
X[,3] <- X[,3]-mean(X[,3])
X <- cbind(1,X)
delta <- as.vector(mydata$delta)
data1=list(Y=Y,delta=delta,X=X)
S1 <- function(data1,theta){
  Y=as.vector(data1$Y)
  X=data1$X
  delta=data1$delta
  z <- 0
  for(i in 1:length(Y)){
    z[i] <- (Y[i]-t(X[i,])%*%theta[-1])/theta[1]
  }
  c1=-r/theta[1]-1/theta[1]*sum(z*delta)+1/theta[1]*sum(z*exp(z))
  c2 <- 0
  for(l in 1:ncol(X)){
    c2[l] <- -1/theta[1]*sum(X[,l]*delta)+1/theta[1]*sum(exp(z)*X[,l])
  }
  c(c1,c2)
}

H1 <- function(data1,theta){
  Y=data1$Y
  X=data1$X
  delta=data1$delta
  z <- 0
  for(i in 1:length(Y)){
    z[i] <- (Y[i]-t(X[i,])%*%theta[-1])/theta[1]
  }
  H12 <- 0
  for(l in 1:ncol(X)){
    H12[l] <- sum(1/theta[1]^2*delta*X[,l])-sum(1/theta[1]^2*exp(z)*X[,l])-sum(1/theta[1]^2*z*exp(z)*X[,l])
  }
  H11 <- r/theta[1]^2+2/theta[1]^2*sum(z*delta)-2/theta[1]^2*sum(z*exp(z))-1/theta[1]^2*sum(z^2*exp(z))
  H22 <- matrix(0,ncol(X),ncol(X))
  for(l in 1:ncol(X)){
    for(s in 1:ncol(X)){
      H22[l,s] <- -1/theta[1]^2*sum(X[,l]*X[,s]*exp(z))
    }
  }
  rbind(c(H11,H12),cbind((H12),H22))
}

#############Newton
eps <- 10
loop <- 0
Loop <- 100
f_Newton <- function(eps, loop, Loop){
  theta.new <- c(0.9,5,0.05,0.01,0.005,0.5,-0.1,-0.8,-0.3)
  theta.old <- theta.new
  while(loop<Loop&eps>1e-15){
    loop <- loop+1
    theta.old <- theta.new
    S <-S1(data1,theta.old)
    H <-H1(data1,theta.old)
    theta.new <- as.vector(theta.old)-solve(H)%*%S
    theta.new <- as.vector(theta.new)
    eps <- sum((theta.old-theta.new)^2)
    eps <- sqrt(eps)
  }
  return(round(theta.new[-2], 3))
}

f_Newton(10, 0, 100)