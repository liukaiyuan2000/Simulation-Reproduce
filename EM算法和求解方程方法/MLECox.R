
CoxSn <- function(data2,theta){
  #this code is used to compute the scores of the cox's model for parametric component
  X <- data2$X
  Y <- data2$Y
  delta <- data2$delta
  R1=data2$R1
  m <- sum(delta)
  n <- dim(X)[1]
  p <- dim(X)[2]
  S=0
  for(i in 1:m){
    S=S+drop(X[i,])-drop(t(X[R1[[i]],])%*%exp(X[R1[[i]],]%*%theta)/sum(exp(X[R1[[i]],]%*%theta)))
    
  }
  S/n
}  
CoxHn <- function(data2,theta){
  #this code is used to compute the Hessian matrix of the cox's model for parametric component
  X <- data2$X
  Y <- data2$Y
  delta <- data2$delta
  R1=data2$R1
  m <- sum(delta)
  n <- dim(X)[1]
  p <- dim(X)[2]
  S=matrix(0,p,p)
  for(i in 1:m){
    index <- drop(X[R1[[i]],]%*%theta)
    temp <-  drop(exp(index))
    Temp <- t(X[R1[[i]],])%*%diag(temp) 
    S=S+ t(X[R1[[i]],])%*%diag(temp)%*%X[R1[[i]],]/sum(temp)-Temp%*%t(Temp)/(sum(temp))^2
    
  }
  -S/n
} 

MLECox <- function(data1){
  Xf <- data1$X
  Yf <- data1$Y
  deltaf<- data1$delta
  p <- dim(X)[2]
  n <- dim(X)[1]
  T0 <- Yf[deltaf==1]
  Xc <- Xf[delta==1,]
  Yc <- Yf[delta==1]
  
  ro <- order(T0)
  T0 <- sort(T0)
  m <- length(T0)
  X <- rbind(Xc[ro,],Xf[deltaf==0,])
  Y <- c(Yc[ro],Yf[deltaf==0])
  delta <- rep(0,n)
  delta[1:m] <- 1
  R1 <- vector("list",m)
  for(i in 1:m)
    R1[[i]] <- which(Y>=Y[i])
  data2=list(Y=Y,X=X,delta=delta,R1=R1)
  
  beta_new <- rep(0,p)
  beta_old <- rep(0,p)
  
  loop <- 0
  Loop <- 5000
  er <- 1
  while(loop<Loop&er>1e-5){
    
    loop <- loop+1
    beta_old <- beta_new
    S=CoxSn(data2,beta_old)
    H=CoxHn(data2,beta_old)
    
    beta_new <- beta_old-solve(H)%*%S
    er <- sqrt(sum((beta_new-beta_old)^2))
    #cat(beta_new,"\n")
  }
  
  h <- rep(0,m)
  for(i in 1:m){
    h[i] <- 1/sum(exp(drop(X[R1[[i]],]%*%beta_new)))
  }
  
  
  out=list(beta.hat=round(beta_new,3),h=round(h,3)
           ,Hn=round(H,3),er=round(er,3),loop=loop)
  out
  
}

Yt=c(6,6,6,6,7,9,10,10,11,13,16,17,19,20,22,23,25,32,32,34,35)
Yc=c(1,1,2,2,3,4,4,5,5,8,8,8,8,11,11,12,12,15,17,22,23)
delta <- c(0,1,1,1,1,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,rep(1,length(Yc)))
Y <- c(Yt,Yc)
X <- c(rep(1,length(Yt)),rep(0,length(Yc)))
X <- cbind(1,X)
data1=list(Y=Y,delta=delta,X=X)
MLECox(data1)
