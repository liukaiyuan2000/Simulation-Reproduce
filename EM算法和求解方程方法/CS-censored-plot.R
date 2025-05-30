setwd("C:/Users/lenovo/Desktop/EM算法和求解方程方法")
 
DATA<-read.table("censor.txt",header = T)
Y <- DATA[,1]
C <- DATA[,2]
X <- DATA[,3]
delta <- DATA[,4]
lambda.new=1/mean(X)
lambda.new1=1/mean(X)
S <- matrix(lambda.new1,50,2)
n=30
for(sim in 2:50){
  lambda.new<-n/(sum(X)+sum(1-delta)/lambda.new)#普通EM更新
  S[sim,1]<-lambda.new
  m=5^(1+floor(sim/10))
  K=0
  for(j in 1:m){
    Z.new=rexp(sum(1-delta),lambda.new1)
    K=K+sum(Z.new)+sum(X)
    
  }
  
  lambda.new1=n*m/K#MCEM更新
  S[sim,2]<-lambda.new1
}

print(S)
plot(S[,1],xlab="index_t",ylab=expression(paste(lambda,"(t)")),type="o",col=4, ylim=c(0.2,0.5))
points(S[,2],type="o",pch=3,col=2,ylim=c(0.2,0.5))
legend(35,0.4,c("EM","MCEM"), col=c(4,2),pch=c(1,3))