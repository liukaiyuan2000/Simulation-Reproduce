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


#############Newton

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
}
print(round(theta.new,3))

cov1=solve(-H)
theta.sd <- sqrt(diag(cov1))
Tn <- theta.new/theta.sd
pvalue <-2*(1-pnorm(Tn))
round(cbind(theta.new,theta.sd,pvalue),3)

#############Gauss_Seidel

theta.hat = c(1,1,1)
itr = 100
theta.values = matrix(0,itr+1,3)
theta.values[1,] = theta.hat
for (j in 1:itr){
  for (i in 1:itr){
    xt = theta.hat[1]-
      ((S1(data1,theta.hat)[1])/(H1(data1,theta.hat)[1,1]))
    theta.values[j+1,1] = xt
    theta.hat[1] = xt
  }
  for (i in 1:itr){
    xt = theta.hat[2]-
      ((S1(data1,theta.hat)[2])/(H1(data1,theta.hat)[2,2]))
    theta.values[j+1,2] = xt
    theta.hat[2] = xt
  }
  for (i in 1:itr){
    xt = theta.hat[3]-
      ((S1(data1,theta.hat)[3])/(H1(data1,theta.hat)[3,3]))
    theta.values[j+1,3] = xt
    theta.hat[3] = xt
  }
}
print(round(theta.hat,3))