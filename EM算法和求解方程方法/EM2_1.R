library("survival")
mydata <- read.csv("可靠性2.1.csv",header=T)
n <- nrow(mydata)
delta <- as.vector(mydata$delta)
mydata1 <- mydata[order(-delta),]
delta <- as.vector(mydata1$delta)
r <- length(delta)-length(delta[delta==0])
st <- mydata1$t
Y <- log(mydata1$t)
X <- as.matrix(mydata1[,2:8])
p=ncol(X)
#################
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
    at = theta.hat[1]-
      ((Q1(theta.hat[1],a_k,
           theta.hat[-1],b_k)[1])
       /(Q2(theta.hat[1],a_k,
            theta.hat[-1],b_k)[1,1]))
    
    theta.hat[1] = at
    bt = theta.hat[2]-
      ((Q1(theta.hat[1],a_k,
           theta.hat[-1],b_k)[2])
       /(Q2(theta.hat[1],a_k,
            theta.hat[-1],b_k)[2,2]))
    
    theta.hat[2] = bt
    ct = theta.hat[3]-
      ((Q1(theta.hat[1],a_k,
           theta.hat[-1],b_k)[3])
       /(Q2(theta.hat[1],a_k,
            theta.hat[-1],b_k)[3,3]))
    
    theta.hat[3] = ct
    
    dt = theta.hat[4]-
      ((Q1(theta.hat[1],a_k,
           theta.hat[-1],b_k)[4])
       /(Q2(theta.hat[1],a_k,
            theta.hat[-1],b_k)[4,4]))
    
    theta.hat[4] = dt
    et = theta.hat[5]-
      ((Q1(theta.hat[1],a_k,
           theta.hat[-1],b_k)[5])
       /(Q2(theta.hat[1],a_k,
            theta.hat[-1],b_k)[5,5]))
    
    theta.hat[5] = et
    ft = theta.hat[6]-
      ((Q1(theta.hat[1],a_k,
           theta.hat[-1],b_k)[6])
       /(Q2(theta.hat[1],a_k,
            theta.hat[-1],b_k)[6,6]))
    
    theta.hat[6] = ft
    
    gt = theta.hat[7]-
      ((Q1(theta.hat[1],a_k,
           theta.hat[-1],b_k)[7])
       /(Q2(theta.hat[1],a_k,
            theta.hat[-1],b_k)[7,7]))
    
    theta.hat[7] = gt
    ht = theta.hat[8]-
      ((Q1(theta.hat[1],a_k,
           theta.hat[-1],b_k)[8])
       /(Q2(theta.hat[1],a_k,
            theta.hat[-1],b_k)[8,8]))
    
    theta.hat[8] = ht
  }
  return(theta.hat)
}


#################################
a_k=0.01
b_k=rep(0.01,p)
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


round(Q1(a_k,a_k,b_k,b_k),3)
round(c(a_k,b_k),4)