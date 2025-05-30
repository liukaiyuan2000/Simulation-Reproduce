library("survival")
library(MASS)
n=50
beta0 <- c(-1,2)
p <- length(beta0)
X=matrix(c(rnorm(n),rnorm(n)),n,p)
U <- runif(n)
alpha0 <- 1
sT <- (-exp(-X%*%beta0)*log(U))^(1/alpha0)
c <- runif(n,0,3)
st=pmin(sT,c)
delta <- 1*(sT<=c)
mean(delta)
r=sum(delta)
mydata=cbind(st,X,delta)
mydata1 <- mydata[order(-delta),]
st <- as.vector(mydata1[,1])
X <- mydata1[,2:(p+1)]
delta <- as.vector(mydata1[,p+2])
###################
res.cox <- coxph(Surv(st,delta)~X)
res.cox$coef


#######log(Ti)的期望
E_lnti <- function(a_k,b_k){
  E_lnti <- rep(0,n)
  E_lnti[1:r] <- log(st[1:r])
  for(i in (r+1):n){
    u <- as.vector(exp(t(X[i,])%*%b_k))
    f1 <- function(x) exp(-x)*log(x/u)
    E_lnti[i] <- 1/a_k*exp(st[i]^a_k*u)*integrate(f1,st[i]^a_k*u,Inf)$value
  }
  #E_lnti <- round(E_lnti,3)
  return(E_lnti)
}


#############Ti^α的期望
E_tia <- function(a,a_k,b_k){
  u <- as.vector(exp(X[(r+1):n,]%*%b_k))
  E_tia <- rep(0,n)
  E_tia[1:r] <- st[1:r]^a
  E_tia[(r+1):n] <- exp(st[(r+1):n]^a_k*u)*as.vector(exp(-a/a_k*X[(r+1):n,]%*%b_k))*
    gamma(a/a_k+1)*(1-pgamma(st[(r+1):n]^a_k*u,a/a_k+1))
  return(E_tia)
}


#############Ti^α的期望对α的一阶导
E_tia1 <- function(a,a_k,b_k){
  f2 <- function(t) t^(a/a_k)*log(t)/exp(t)
  E_tia1 <- rep(0,n)
  E_tia1[1:r] <- st[1:r]^a*log(st[1:r])
  for(i in (r+1):n){
    u <- as.vector(exp(t(X[i,])%*%b_k))
    E_tia1[i] <- -1/a_k*as.vector(t(X[i,])%*%b_k)*E_tia(a,a_k,b_k)[i]+
      exp(st[i]^a_k*u)*as.vector(exp(-a/a_k*t(X[i,])%*%b_k))*
      integrate(f2,st[i]^a_k*u,Inf)$value/a_k
  }
  #E_tia1 <- round(E_tia1,3)
  return(E_tia1)
}


#############Ti^α的期望对α的二阶导
E_tia2 <- function(a,a_k,b_k){
  f2 <- function(t) t^(a/a_k)*log(t)/exp(t)
  f3 <- function(t) t^(a/a_k)*(log(t))^2/exp(t)
  E_tia2 <- rep(0,n)
  E_tia2[1:r] <- st[1:r]^a*(log(st[1:r]))^2
  for(i in (r+1):n){
    u <- as.vector(exp(t(X[i,])%*%b_k))
    E_tia2[i] <- -1/a_k*as.vector(t(X[i,])%*%b_k)*
      E_tia1(a,a_k,b_k)[i]-1/a_k^2*exp(st[i]^a_k*u)*
      as.vector(t(X[i,])%*%b_k)*as.vector(exp(-a/a_k*t(X[i,])%*%b_k))*
      integrate(f2,st[i]^a_k*u,Inf)$value+exp(st[i]^a_k*u)*
      as.vector(exp(-a/a_k*t(X[i,])%*%b_k))*
      integrate(f3,st[i]^a_k*u,Inf)$value/a_k^2
  }
  #E_tia2 <- round(E_tia2,3)
  return(E_tia2)
}


#############似然函数Ln(l(θ))的期望
Q <- function(a,a_k,b,b_k){
  Q <- 0
  for(i in 1:n){
    Q=Q+as.vector(t(X[i,])%*%b)-as.vector(exp(t(X[i,])%*%b))*
      E_tia(a,a_k,b_k)[i]+log(a)+(a-1)*E_lnti(a_k,b_k)[i]
  }
  return(Q)
}


#############似然函数Ln(l(θ))的期望对θ的一阶导
Q1 <- function(a,a_k,b,b_k){
  Q1_a <- sum(-E_tia1(a,a_k,b_k)*as.vector(exp(X%*%b))+1/a+E_lnti(a_k,b_k))
  Q1_b <- apply((X-X*as.vector(exp(X%*%b))*E_tia(a,a_k,b_k)),2,sum)
  Q1=c(Q1_a,Q1_b)
  return(Q1)
}


#############似然函数Ln(l(θ))的期望对θ的二阶导
Q2 <- function(a,a_k,b,b_k){
  Q2_a <- sum(-E_tia2(a,a_k,b_k)*as.vector(exp(X%*%b))-1/a^2)
  Q2_b <- -t(X*(as.vector(exp(X%*%b))*E_tia(a,a_k,b_k)))%*%X
  Q2_ab <- apply((-X*as.vector(exp(X%*%b))*E_tia1(a,a_k,b_k)),2,sum)
  Q2 <- rbind(c(Q2_a,Q2_ab),cbind(Q2_ab,Q2_b))
  return(Q2)
}


###########################Gauss-Seidal迭代
guass_seidal <- function(theta.hat){
  loop=0
  while(sqrt(sum(Q1(theta.hat[1],a_k,theta.hat[-1],b_k)^2))>1e-4&loop<100){
    loop=loop+1
    cat(loop,"\r")
    xt = theta.hat[1]-
      ((Q1(theta.hat[1],a_k,
           theta.hat[-1],b_k)[1])
       /(Q2(theta.hat[1],a_k,
            theta.hat[-1],b_k)[1,1]))
    
    theta.hat[1] = xt
    yt = theta.hat[2]-
      ((Q1(theta.hat[1],a_k,
           theta.hat[-1],b_k)[2])
       /(Q2(theta.hat[1],a_k,
            theta.hat[-1],b_k)[2,2]))
    
    theta.hat[2] = yt
    zt = theta.hat[3]-
      ((Q1(theta.hat[1],a_k,
           theta.hat[-1],b_k)[3])
       /(Q2(theta.hat[1],a_k,
            theta.hat[-1],b_k)[3,3]))
    
    theta.hat[3] = zt
  }
  return(theta.hat)
}


############################单次求解
a_k=0.8
b_k=c(-0.8,1.8)
theta_k <- matrix(0,50,p+1)
theta_k[1,] <- c(a_k,b_k)
k=1
while(sqrt(sum(Q1(a_k,a_k,b_k,b_k)^2))>1e-4&k<50){
  k=k+1
  cat(k,"\r\n")
  theta.hat <- theta_k[(k-1),]
  theta_k[k,] <- guass_seidal(theta.hat)
  a_k <- theta_k[k,1]
  b_k <- theta_k[k,-1]
}
theta_k


###############Bootstrap--痛点在于复杂度高导致循环时间过长
B=20
myroot1 <- matrix(0,B,p+1)
i=1
while(i<=B){
  n=50
  p <- length(beta0)
  X=matrix(c(rnorm(n),rnorm(n)),n,p)
  U <- runif(n)
  alpha0 <- 1
  sT <- (-exp(-X%*%beta0)*log(U))^(1/alpha0)
  c <- runif(n,0,6)#--设置删失率(3-40%,6-30%,12-20%)
  st=pmin(sT,c)
  delta <- 1*(sT<=c)
  if(mean(delta)>0.69&mean(delta)<0.71){#--设置删失率
    r=sum(delta)
    mydata=cbind(st,X,delta)
    mydata1 <- mydata[order(-delta),]
    st <- as.vector(mydata1[,1])
    X <- mydata1[,2:(p+1)]
    delta <- as.vector(mydata1[,p+2])
    #初值设置
    a_k=0.8
    b_k=c(-0.8,1.8)
    theta_k <- matrix(0,50,p+1)
    theta_k[1,] <- c(a_k,b_k)
    k=1
    while(sqrt(sum(Q1(a_k,a_k,b_k,b_k)^2))>1e-4&k<50){
      k=k+1
      cat(i," ",k,"\r\n")
      theta.hat <- theta_k[(k-1),]
      theta_k[k,] <- guass_seidal(theta.hat)
      a_k <- theta_k[k,1]
      b_k <- theta_k[k,-1]
    }
    myroot1[i,] <- c(a_k,b_k)
    i=i+1
  }
  else i=i
}
BEA_1 <- sum(myroot1[,1]-alpha0)/B#BE-alpha
BEB1_1 <- sum(myroot1[,2]-beta0[1])/B#BE-beta1
BEB2_1 <- sum(myroot1[,3]-beta0[2])/B#BE-beta2
MSEA_1 <- sum((myroot1[,1]-alpha0)^2)/B#MSE-alpha
MSEB1_1 <- sum((myroot1[,2]-beta0[1])^2)/B#MSE-beta1
MSEB2_1 <- sum((myroot1[,3]-beta0[2])^2)/B#MSE-beta2
round(c(BEA_1,BEB1_1,BEB2_1),3)
round(c(MSEA_1,MSEB1_1,MSEB2_1),3)
round(apply(myroot1,2,mean),3)


####################################
myroot2 <- matrix(0,B,p+1)
i=1
while(i<=B){
  n=100
  p <- length(beta0)
  X=matrix(c(rnorm(n),rnorm(n)),n,p)
  U <- runif(n)
  alpha0 <- 1
  sT <- (-exp(-X%*%beta0)*log(U))^(1/alpha0)
  c <- runif(n,0,3)#--设置删失率(3-40%,6-30%,12-20%)
  st=pmin(sT,c)
  delta <- 1*(sT<=c)
  if(mean(delta)>0.58&mean(delta)<0.62){#--设置删失率
    r=sum(delta)
    mydata=cbind(st,X,delta)
    mydata1 <- mydata[order(-delta),]
    st <- as.vector(mydata1[,1])
    X <- mydata1[,2:(p+1)]
    delta <- as.vector(mydata1[,p+2])
    #初值设置
    a_k=0.8
    b_k=c(-0.8,1.8)
    theta_k <- matrix(0,50,p+1)
    theta_k[1,] <- c(a_k,b_k)
    k=1
    while(sqrt(sum(Q1(a_k,a_k,b_k,b_k)^2))>1e-4&k<50){
      k=k+1
      cat(i," ",k,"\r\n")
      theta.hat <- theta_k[(k-1),]
      theta_k[k,] <- guass_seidal(theta.hat)
      a_k <- theta_k[k,1]
      b_k <- theta_k[k,-1]
    }
    myroot2[i,] <- c(a_k,b_k)
    i=i+1
  }
  else i=i
}
BEA_2 <- sum(abs(myroot2[,1]-alpha0))/B#BE-alpha
BEB1_2 <- sum(abs(myroot2[,2]-beta0[1]))/B#BE-beta1
BEB2_2 <- sum(abs(myroot2[,3]-beta0[2]))/B#BE-beta2
MSEA_2 <- sum((myroot2[,1]-alpha0)^2)/B#MSE-alpha
MSEB1_2 <- sum((myroot2[,2]-beta0[1])^2)/B#MSE-beta1
MSEB2_2 <- sum((myroot2[,3]-beta0[2])^2)/B#MSE-beta2
round(c(BEA_2,BEB1_2,BEB2_2),3)
round(c(MSEA_2,MSEB1_2,MSEB2_2),3)
round(apply(myroot2,2,mean),3)


##################################
myroot3 <- matrix(0,B,p+1)
i=1
while(i<=B){
  n=200
  p <- length(beta0)
  X=matrix(c(rnorm(n),rnorm(n)),n,p)
  U <- runif(n)
  alpha0 <- 1
  sT <- (-exp(-X%*%beta0)*log(U))^(1/alpha0)
  c <- runif(n,0,6)#--设置删失率(3-40%,6-30%,12-20%)
  st=pmin(sT,c)
  delta <- 1*(sT<=c)
  if(mean(delta)>0.68&mean(delta)<0.72){#--设置删失率
    r=sum(delta)
    mydata=cbind(st,X,delta)
    mydata1 <- mydata[order(-delta),]
    st <- as.vector(mydata1[,1])
    X <- mydata1[,2:(p+1)]
    delta <- as.vector(mydata1[,p+2])
    #初值设置
    a_k=0.8
    b_k=c(-0.8,1.8)
    theta_k <- matrix(0,50,p+1)
    theta_k[1,] <- c(a_k,b_k)
    k=1
    while(sqrt(sum(Q1(a_k,a_k,b_k,b_k)^2))>1e-4&k<50){
      k=k+1
      cat(i," ",k,"\r\n")
      theta.hat <- theta_k[(k-1),]
      theta_k[k,] <- guass_seidal(theta.hat)
      a_k <- theta_k[k,1]
      b_k <- theta_k[k,-1]
    }
    myroot3[i,] <- c(a_k,b_k)
    i=i+1
  }
  else i=i
}
BEA_3 <- sum(abs(myroot3[,1]-alpha0))/B#BE-alpha
BEB1_3 <- sum(abs(myroot3[,2]-beta0[1]))/B#BE-beta1
BEB2_3 <- sum(abs(myroot3[,3]-beta0[2]))/B#BE-beta2
MSEA_3 <- sum((myroot3[,1]-alpha0)^2)/B#MSE-alpha
MSEB1_3 <- sum((myroot3[,2]-beta0[1])^2)/B#MSE-beta1
MSEB2_3 <- sum((myroot3[,3]-beta0[2])^2)/B#MSE-beta2
round(c(BEA_3,BEB1_3,BEB2_3),3)
round(c(MSEA_3,MSEB1_3,MSEB2_3),3)
round(apply(myroot3,2,mean),3)


################################
myroot4 <- matrix(0,B,p+1)
i=1
while(i<=B){
  n=400
  p <- length(beta0)
  X=matrix(c(rnorm(n),rnorm(n)),n,p)
  U <- runif(n)
  alpha0 <- 1
  sT <- (-exp(-X%*%beta0)*log(U))^(1/alpha0)
  c <- runif(n,0,3)#--设置删失率(3-40%,6-30%,12-20%)
  st=pmin(sT,c)
  delta <- 1*(sT<=c)
  if(mean(delta)>0.58&mean(delta)<0.62){#--设置删失率
    r=sum(delta)
    mydata=cbind(st,X,delta)
    mydata1 <- mydata[order(-delta),]
    st <- as.vector(mydata1[,1])
    X <- mydata1[,2:(p+1)]
    delta <- as.vector(mydata1[,p+2])
    #初值设置
    a_k=0.8
    b_k=c(-0.8,1.8)
    theta_k <- matrix(0,50,p+1)
    theta_k[1,] <- c(a_k,b_k)
    k=1
    while(sqrt(sum(Q1(a_k,a_k,b_k,b_k)^2))>1e-4&k<50){
      k=k+1
      cat(i," ",k,"\r\n")
      theta.hat <- theta_k[(k-1),]
      theta_k[k,] <- guass_seidal(theta.hat)
      a_k <- theta_k[k,1]
      b_k <- theta_k[k,-1]
    }
    myroot4[i,] <- c(a_k,b_k)
    i=i+1
  }
  else i=i
}
BEA_4 <- sum(abs(myroot4[,1]-alpha0))/B#BE-alpha
BEB1_4 <- sum(abs(myroot4[,2]-beta0[1]))/B#BE-beta1
BEB2_4 <- sum(abs(myroot4[,3]-beta0[2]))/B#BE-beta2
MSEA_4 <- sum((myroot4[,1]-alpha0)^2)/B#MSE-alpha
MSEB1_4 <- sum((myroot4[,2]-beta0[1])^2)/B#MSE-beta1
MSEB2_4 <- sum((myroot4[,3]-beta0[2])^2)/B#MSE-beta2
round(c(BEA_4,BEB1_4,BEB2_4),3)
round(c(MSEA_4,MSEB1_4,MSEB2_4),3)
round(apply(myroot4,2,mean),3)