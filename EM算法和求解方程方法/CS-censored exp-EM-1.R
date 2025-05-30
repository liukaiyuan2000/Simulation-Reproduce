rm(list=ls())  
set.seed(1234567)
Sim=2000
n=500
lambda=1/200
Sim=1000
SD=matrix(0,Sim,4)
for(sim in 1:Sim){
  
  Y=rexp(n,lambda)
  lambda.gold=mean(Y)
  #C=runif(n,150,180)
  C=rexp(n,lambda/2)
  delta=1*(Y<C)
  mean(delta)
  X=pmin(Y,C)
  lambda.hat=1/(sum(X)/sum(delta))
  LB=NULL
  B <- 1000
  for(b in 1:B){
    index=sample.int(n,n,replace = T)
    X.new=X[index]
    delta.new=delta[index]
    lambda.b=1/mean(X.new)
    for(it in 1:20)
      lambda.b=n/(sum(X.new)+(n-sum(delta.new))/ lambda.b)
    LB=c(LB,lambda.b)
    
  }
  sd1=sqrt(1/(sum(delta)/lambda.hat^2))   
  sd2=sd(LB)
  Low= lambda.hat-1.96*sd1
  Up= lambda.hat+1.96*sd1
  ecp1=1*(lambda<Up)*(lambda>Low)#落在置信区间的概率
  Low= lambda.hat-1.96*sd2
  Up= lambda.hat+1.96*sd2
  ecp2=1*(lambda<Up)*(lambda>Low)#落在置信区间的概率
  SD[sim,]=c(sd1,ecp1,sd2,ecp2)
  
}
 out=round(rbind(apply(SD,2,mean),apply(SD,2,sd))*100,2)
 print(out)
