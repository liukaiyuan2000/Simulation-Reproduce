
############################
#computational statistics 
##Exercise 2.3##
rm(list=ls())
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
  H12 <- apply( diag(as.vector(Y^theta[1]*log(Y)*exp(-index)))%*%X ,2,sum)
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
  ###cat(c(eps,theta.new),"\n")
}
out1=list(theta.new,loop)

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
  H12 <- apply( diag(as.vector(Y^theta[1]*log(Y)*exp(-index)))%*%X ,2,sum)
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
out1=list(theta.new=theta.new,loop=loop)

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


 
##########################################
#########################################################################
### Exercise 2.3 STEEPEST ASCENT
#########################################################################
Ln <- function(data1,theta){
  #theta <- theta.new
  Y=data1$Y
  X=data1$X
  delta=data1$delta
  index <- as.vector(X%*%theta[-1])
  Lammda <- as.vector(Y^theta[1])
  lammda <- as.vector(theta[1]*Y^(theta[1]-1))
  logSt <- -Lammda*exp(-index)
  htx <- lammda*exp(-index)
  sum(delta*log(htx)+logSt)
  
}

x = c(1,3,1)
M = -diag(p)
itr = 500
alpha.default = 0.01
alpha = alpha.default
x.values = matrix(0,itr+1,p)
x.values[1,] = x



## MAIN
for (i in 1:itr){
  hessian.inv = solve(M)
  ht <-  -alpha*hessian.inv%*%S1(data1,x)
  xt = x +ht
  # REDUCE ALPHA UNTIL A CORRECT STEP IS REACHED
  while(Ln(data1,xt) < Ln(data1,x)){
    alpha = alpha/2
    xt = x - alpha*hessian.inv%*%S1(data1,x)
    
  }
  x.values[i+1,] = x = xt
  alpha = alpha.default
}
 
out2=list(theta.new=x)

##nonlinear Guass-Seidel iteration

NR <- function(data1,theta){
  Y=data1$Y
  X=data1$X
  delta=data1$delta
  index <- as.vector(X%*%theta[-1])
  loop1 <- 0
  Loop1 <- 100
  er1 <- 1
  while(loop1<Loop1&er1>1e-5){
    loop1 <- loop1+1
    theta.old <- theta.new
    SS <- S1(data1,theta.old)
    SS.temp <- SS[1]
    HH <-as.vector(H1(data1,theta.old))
    HH.temp <- HH[1]
    theta.new[1] <- theta.old[1]- SS.temp/HH.temp
    
    er <- abs(theta.new[1]- theta.old[1])
  }
  theta.new
  
}




NRI <- function(data1,theta,i){
  Y=data1$Y
  X=data1$X
  delta=data1$delta
  index <- as.vector(X%*%theta[-1])
  loop1 <- 0
  Loop1 <- 100
  er1 <- 1
  while(loop1<Loop1&er1>1e-5){
    loop1 <- loop1+1
    theta.old <- theta.new
    SS <- S1(data1,theta.old)
    SS.temp <- SS[i]
    HH <- H1(data1,theta.old) 
    HH.temp <- HH[i,i]
    theta.new[i] <- theta.old[i]- SS.temp/HH.temp
    
    er <- abs(theta.new[i]- theta.old[i])
  }
  theta.new
  
}



loop <- 0
Loop <- 100
Y=data1$Y
X=data1$X
delta=data1$delta
p <- ncol(X)+1
theta.new <- runif(p)
theta.old <- runif(p)

 
er <- 1
while(loop<Loop&er>1e-5){
  loop <- loop+1
  theta.old <- theta.new
  for(i in 1:p){
    if(i==1) theta.new  <- NR(data1,theta.new)
    if(i>1) theta.new  <- NRI(data1,theta.new,i)
  }
  eps <- sum((theta.old - theta.new)^2)
  er <- sqrt(eps)
  
}


DN <- function(data1){
  Y=data1$Y
  X=data1$X
  delta=data1$delta
  p <- ncol(X)+1
  n <- nrow(X)
  theta.new <- runif(p)
  theta.old <- runif(p)
  loop1 <- 0
  Loop1 <- 100
  er1 <- 1
  h <- (1e-3)/n
  while(loop1<Loop1&er1>1e-5){
    loop1 <- loop1+1
    theta.old <- theta.new
    S <- S1(data1,theta.old)
    M <- matrix(0,p,p)
    for(i in 1:p){ 
      theta.temp <- theta.old
      theta.temp[i] <- theta.temp[i]+h
      M[i,] <-  (S1(data1,theta.temp)-S)/h
    }
       
    theta.new <- as.vector(theta.old)-solve(M)%*%S
    theta.new <- as.vector( theta.new)
     
    
    er1 <- sum((theta.old - theta.new)^2)
    er1 <- sqrt(er1)
    #cat(c(theta.new=theta.new,loop=loop,er1),"\n")
  }
  out=list(theta.new=theta.new,loop=loop,er1)
  out
}

out3=DN(data1)

out=round(cbind(out1$theta.new,out2$theta.new,theta.new,out3$theta.new),3)

############################
#computational statistics 
##Exercise 2.4##
alpha <- 0.05

f <- function(x) x*exp(-x)+(exp(-x)-1+alpha)*log(exp(-x)-1+alpha)
f1 <- function(x) -x*exp(-x)-exp(-x)*log(exp(-x)-1+alpha)

x <- 0.003
for(i in 1:100){
  x <- x-f(x)/f1(x)
  cat(x,"\n")
  
}
 
y <- -log(-1+alpha+exp(-x))
a=1-x*exp(-x)-exp(-x)
b=1-y*exp(-y)-exp(-y)
b-a

x*exp(-x)-y*exp(-y)

t <- seq(0.01,0.051,length=200)
par(mfrow = c(1, 2))

plot(t,f(t),type="l",col=4)
#par(new=TRUE)
plot(t,f1(t),type="l",col=4)
 
 



rm(list=ls())

Pf <- function(data1,theta){
  #theta=runif(2)
  Y=data1$Y
  t=data1$t
  N0=data1$N0
  c1=N0^2*(1-exp(-theta[2]*t))/(N0+(theta[1]-N0)*exp(-theta[2]*t))/(N0+(theta[1]-N0)*exp(-theta[2]*t))  
  c2=N0*theta[1]*(theta[1]-N0)*exp(-theta[2]*t)*t/(N0+(theta[1]-N0)*exp(-theta[2]*t))/(N0+(theta[1]-N0)*exp(-theta[2]*t))
  cbind(c1,c2)
}
ft <- function(data1,theta){
  Y=data1$Y
  t=data1$t
  N0=data1$N0
  N0*theta[1]/(N0+(theta[1]-N0)*exp(-theta[2]*t))
  
}

setwd("C:/Users/Administrator/Desktop/CompStat/Code")
data1 <- read.table("C:/Users/Administrator/Desktop/CompStat/mydata/flourbeetles.txt", header = TRUE, sep = "", quote = "\"'")

theta.new <- runif(2)
theta.old <- runif(2)
t <- data1[,1]
Y <- data1[,2]
 data1 =list(t=t,Y=Y,N0=Y[1])
loop1 <- 0
Loop1 <- 1000
er1 <- 1
 
while(loop1<Loop1&er1>1e-5){
  loop1 <- loop1+1
  theta.old <- theta.new
  Xt <-Y- ft(data1,theta.old)
  A <- Pf(data1,theta.old)
  theta.new <- theta.old+solve(t(A)%*%A+1e-7)%*%t(A)%*%Xt 
  er1 <- sum((theta.old - theta.new)^2)
  er1 <- sqrt(er1)
  cat(c(theta.new,loop1,er1),"\n")
} 


mean((Y-ft(data1,theta.new))^2)



S1 <- function(data1,theta){
  #theta=runif(2)
  Y=data1$Y
  t=data1$t
  N0=data1$N0
  c1=N0^2*(1-exp(-theta[2]*t))/(N0+(theta[1]-N0)*exp(-theta[2]*t))/(N0+(theta[1]-N0)*exp(-theta[2]*t))  
  c2=N0*theta[1]*(theta[1]-N0)*exp(-theta[2]*t)*t/(N0+(theta[1]-N0)*exp(-theta[2]*t))/(N0+(theta[1]-N0)*exp(-theta[2]*t))
  
  t(Y-ft(data1,theta))%*%cbind(c1,c2)
}

NH <- function(data1,theta){
  #theta=theta.old
  Y=data1$Y
  n <- length(Y)
  t=data1$t
  N0=data1$N0
  c1=N0^2*(1-exp(-theta[2]*t))/(N0+(theta[1]-N0)*exp(-theta[2]*t))/(N0+(theta[1]-N0)*exp(-theta[2]*t))  
  c2=N0*theta[1]*(theta[1]-N0)*exp(-theta[2]*t)/(N0+(theta[1]-N0)*exp(-theta[2]*t))/(N0+(theta[1]-N0)*exp(-theta[2]*t))
  temp1 <- -t(cbind(c1,c2))%*% cbind(c1,c2)
  res=(Y-ft(data1,theta))
  p <- length(theta)
  H <- matrix(0,p,p)
  for(i in 1:n){
    H <- H+res[i]*DL(N0,t[i],theta)
  }
  temp1+H
  
}

DL <- function(N0,t,theta){ 
  
  d1 <- N0+(theta[1]-N0)*exp(-theta[2]*t)
  f11 <- -2*N0^2*(1-exp(-theta[2]*t))*exp(-theta[2]*t)*(d1^(-3))
  f12 <- N0^2*exp(-theta[2]*t)*t*(d1^(-2))+2*N0^2*theta[1]*(1-exp(-theta[2]*t))*d1^(-3)*exp(-theta[2]*t)*t
  f21 <- f12
  f22 <- -N0*theta[1]*(theta[1]-N0)*exp(-theta[2]*t)*t^2*d1^(-2)+2*N0*theta[1]*(theta[1]-N0)*exp(-theta[2]*t)*(theta[1]-N0)*exp(-theta[2]*t)*t^2*d1^(-3)
  matrix(c(f11,f12,f21,f22),2,2)
}

ft <- function(data1,theta){
  Y=data1$Y
  t=data1$t
  N0=data1$N0
  N0*theta[1]/(N0+(theta[1]-N0)*exp(-theta[2]*t))
  
}





theta.new <- runif(2)
theta.old <- runif(2)
 
loop1 <- 0
Loop1 <- 1000
er1 <- 1

while(loop1<Loop1&er1>1e-5){
  loop1 <- loop1+1
  theta.old <- theta.new
   SS <- S1(data1,theta.old)
   HH <- NH(data1,theta.old)
   theta.new <- theta.old-solve(HH)%*%as.vector(SS)
  er1 <- sum((theta.old - theta.new)^2)
  er1 <- sqrt(er1)
  cat(c(theta.new,loop1,er1),"\n")
} 


mean((Y-ft(data1,theta.new))^2)





