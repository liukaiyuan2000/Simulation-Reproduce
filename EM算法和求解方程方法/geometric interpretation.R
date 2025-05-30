
x=seq(1,5,length=1000)
y=log(x)/(1+x)
y1=(1+1/x-log(x))/(1+x)
y1=y1/(1+x)
par(mfrow = c(1, 2))
plot(x,y,type="l",col=2, main = "the plot of g(x)") 
plot(x,y1,type="l",col=3, main = "the first derivative of g(x)")




x=seq(0.5,1.5,length=1000)
g=(1-x)^2/10000
f=x^3*(1-x)^2
par(mfrow = c(1, 2))
plot(x,g,type="l",col=2, main = "the plot of g(x)") 
plot(x,f,type="l",col=3, main = "the plot of f(x)") 


#函数f(x)与g(x)
x=seq(-0.5,0.5,length=1000)
g=dnorm(x,0,100)
f=dnorm(x,0,0.01)
par(mfrow = c(1, 2))
plot(x,g,type="l",col=2, main = "the plot of g(x)") 
plot(x,f,type="l",col=3, main = "the plot of f(x)") 






 #用0.618方法求闭区间[a,b]上单峰函数的最小值 

 ZSOE<-function(f,a0=1e-5,b0=1,eps=1e-8,Loop=100){
    tol = .Machine$double.eps^0.5
    tau=sqrt(5)/2-0.5    
    loop=0
    alpha_r=a0+tau*(b0-a0)
    alpha_l=a0+(1-tau)*(b0-a0)     
    while((b0-a0>eps)&(loop<Loop)){
      loop=loop+1
      alpha_r=a0+tau*(b0-a0)
      alpha_l=a0+(1-tau)*(b0-a0)  
      if(f(alpha_l)<f(alpha_r)) b0=alpha_r else  a0=alpha_l
      alpha=(a0+b0)/2
     
    }   
    list(root=alpha,loop,eps=b0-a0)
  }

  f<-function(x) 2*x^2-x-1
 
  ZSOE(f,a0=0.1,b0=0.5)
  x=seq(-10,10,length=1000) 
  plot(x,sapply(x,f),type="l",col=3)




#二分法求根
dichotomy<-function(f,a,b,eps=1e-5){
  if(f(a)*f(b)>0)
    list(fail="find root is fail!")
  else{
    repeat{
      if(abs(b-a)<eps) break;
      x<-(a+b)/2
      if(f(a)*f(x)<0) b<-x else a<-x
    }
    list(root=(a+b)/2,fun=f(x))
  }
}
f<-function(x) x^3-x-1
#f<-function(x) (1+1/x-log(x))/(1+x)
dichotomy(f,1,2)
dichotomy(f,0,1)



#二分法求根
dichotomy<-function(f,a,b,eps=1e-5,Loop=100){
  loop=0
  x_old=0
  x_new=1
  if(f(a)*f(b)>0)
    list(fail="find root is fail!")
  else{
    while((abs(x_new-x_old)>1e-29)&loop<Loop){
      loop=loop+1
      x_old=x_new
      x_new=(a+b)/2
      if(f(a)*f(b)>0) break;     
      if(f(a)*f(x_new)<0) b<-x_new else a<-x_new
    }   
    list(root=(a+b)/2,fun=f(x_new),esp=b-a,loop)
  }
}

f<-function(x) x^3-x-1
dichotomy(f,-10,20,eps=1e-5)



 

  #Newton—Raphson
  options(digits=8)
  Newton<-function(f,f1,x0,eps=1e-5,Loop=100){
    tol = .Machine$double.eps^0.5
    loop=0
    x_old=0
    x_new=x0
    while((abs(x_new-x_old)>1e-5)&loop<Loop&abs(f(x_new))>1e-5){
      loop=loop+1
      x_old=x_new
      L1=f(x_old)
      L2=f1(x_old)
      x_new=x_old-L1/(L2+tol)
    }   
    list(root=x_new,esp=abs(x_new-x_old),loop,f0=f(x_new),hessian=f1(x_new))
  }
  f<-function(x) x^3-x-1
  f1<-function(x) 3*x^2-1
  out=Newton(f,f1,2)
  out$root
############################牛顿算法的解释


f<-function(x) 0.25*(x-1)^4+exp(x)
f1<-function(x) (x-1)^3+exp(x)
f2<-function(x) 3*(x-1)^2+exp(x)
q<-function(x,x0) f(x0)+f1(x0)*(x-x0)+0.5*f2(x0)*(x-x0)^2
t=seq(-1,3.5,length=2000)
x0=2.5
t1=seq(1,3,length=200)  
plot(t,f(t),type='l',ylim=c(0,45),xlab="x",ylab="f(x)",xlim=c(-1,4),bty="L")
par(new=TRUE)
plot(t1, q(t1,x0),type='l',ylim=c(0,45),xlab="x",ylab="f(x)",col=2,xlim=c(-1,4))


x1=x0-f1(x0)/f2(x0) 
t2=seq(x1-1.2,x1+0.5,length=200)
par(new=TRUE)
plot(t2, q(t2,x1),type='l',ylim=c(0,45),xlab="x",ylab="f(x)",col=3,xlim=c(-1,4))

x2=x1-f1(x1)/f2(x1)
par(new=TRUE)
t3=seq(x2-1,x2+1,length=200) 
plot(t3, q(t3,x3),type='l',ylim=c(0,45),xlab="x",ylab="f(x)",col=5,xlim=c(-1,4))

 
x3=x2-f1(x2)/f2(x2)
lines(rep(x0,200),seq(-0.10,f(x0),length=200),lty=2,col=4,lwd=1)
lines(rep(x1,200),seq(-0.10,f(x1),length=200),lty=2,col=4,lwd=1)
lines(rep(x2,200),seq(-0.10,f(x2),length=200),lty=2,col=4,lwd=1)
lines(rep(x3,200),seq(-0.10,f(x3),length=200),lty=2,col=4,lwd=1)

text(x0,0.05,expression(x[0]==2.5))
text(x1,-0.05,expression(x[1]),col=2)
text(x2,-0.05,expression(x[2]),col=3)
text(x3,-0.05,expression(x[3]),col=3)

###########################


rm(list = ls())
f<-function(x) 1.95-exp(-2/x)-2*exp(-x^4)
f1<-function(x0) -2*x0^(-2)*exp(-2/x0)+8*x0^3*exp(-x0^4)
 out=Newton(f,f1,1)
  out$root


x=seq(0.2,2,by=0.01)

plot(x,f(x),type='l',lwd=2)
abline(h=0,v=0)
#########
x0=1
gx0=f(x0)
dgx0=f1(x0)
x1=x0-gx0/dgx0
x=seq(x1,x0,by=0.01)
y1=-dgx0*x1+dgx0*x
lines(x,y1,lty=2,col=2,lwd=1)

xy0=seq(0,gx0,by=0.01)
xx0=rep(x0,length(xy0))
lines(xx0,xy0,col=2,lwd=0.15)
#########

gx1=f(x1)
dgx1=f1(x1)
x2=x1-gx1/dgx1
x=seq(x2,x1,by=0.01)
y2=-dgx1*x2+dgx1*x
lines(x,y2,lty=2,col=3,lwd=1)

xy1=seq(0,gx1,by=0.01)
xx1=rep(x1,length(xy1))
lines(xx1,xy1,col=3,lwd=1)

text(x0,-0.05,expression(x[0]==1))
text(x1,-0.05,expression(x[1]),col=2)
text(x2,-0.05,expression(x[2]),col=3)

text(1.7,1.6,expression(g(x)))

#########
gx2=f(x2)
dgx2=f1(x2)
x3=x2-gx2/dgx2
x=seq(x3,x2,by=0.01)
y3=-dgx2*x3+dgx2*x
lines(x,y3,lty=2,col=4,lwd=1)

xy2=seq(0,gx2,by=0.01)
xx2=rep(x2,length(xy2))
lines(xx2,xy2,col=4,lwd=1)
text(x3-0.02,-0.05,expression(x[3]),col=4)





 secant1<-function(f,x0,eps=1e-5,Loop=100){
    tol = .Machine$double.eps^0.5
    loop=0
    x_old=0
    x_new=x0
    while((abs(x_new-x_old)>1e-5)&loop<Loop&abs(f(x_new))>1e-5){
      loop=loop+1
       
      x_old=x_new
      L1=f(x_old)
      L2=(f(x_old)-f(x_old-eps))/eps
      x_new=x_old-L1/(L2+tol)
    }   
    list(root=x_new,esp=abs(x_new-x_old),loop,f0=f(x_new))
  }


###########################切线法的几何解释
 f<-function(t) log(t)/(1+t)
 f1<-function(t) (1+1/t-log(t))/(1+t) 
 x=seq(1.5,5,length=1000)
 
 X=matrix(,10,3)
 X[1,2:3]=c(1.5,2)
 for(i in 2:10){   
   temp=X[i-1,3]-f1(X[i-1,3])*(X[i-1,3]-X[i-1,2])/(f1(X[i-1,3])-f1(X[i-1,2]))
   X[i,]=c(X[i-1,2],X[i-1,3],temp)
 } 

 x1 <- c(1.5, 5)
 y1 <- c(0,0)
 plot(x,f1(x),type="l",col=1,ylab=expression(f(x)))
 lines(x1, y1, lty=3) 
 for(i in 2:6){
  points(X[i, 1], f1(X[i, 1]), pch=20,col=i)
 }

 for(i in 2:5){
   xt <- X[i, 3]+0.1
   yt <- f1(X[i, 1])+(f1(X[i+1, 1])-f1(X[i, 1]))/(X[i+1, 1]-X[i, 1])*(xt-X[i, 1])
   xtt <- X[i, 1]-0.2
   ytt <- f1(X[i, 1])+(f1(X[i+1, 1])-f1(X[i, 1]))/(X[i+1, 1]-X[i, 1])*(xtt-X[i, 1])
   x1 <- c(xtt, xt)
   y1 <- c(ytt, yt)
   lines(x1, y1, lty=3, lwd=1, col=i)
 }

 for(i in 2:4){
   x1 <- c(X[i, 3], X[i, 3])
   y1 <- c(f1(X[i, 3]), 0)
   lines(x1, y1, lty=3, lwd=1.5)
 }

 text(X[2, 1], -0.02, expression(x^{(0)}))
 text(X[2, 2], -0.02, expression(x^{(1)}))
 text(X[2, 3], -0.02, expression(x^{(2)}))
 text(X[3, 3], -0.02, expression(x^{(3)}))
 text(X[4, 3], -0.02, expression(x^{(4)}))

#################################
 #########用压缩映射原理求非线性方程的零点
 #求方程x^2-3=0的根
 #g1(x)=x+x^2-3  #不满足理论要求
 #g2(x)=3/x      #不满足理论要求
 #g3(x)=x-(x^2-3)/4  #满足理论要求
 #g4(x)=(x+3/x)/2 #满足理论要求，最优的
 #g5<-function(x) x+(x^2-3)/100

 g1<-function(x) x+x^2-3
 g2<-function(x) 3/x
 g3<-function(x) x-(x^2-3)/4
 g4<-function(x) (x+3/x)/2
 g5<-function(x) x+(x^2-3)/100
 Loop=10
 X=matrix(1,Loop,5)
 for(i in 2:Loop){
   temp=X[i-1,]
   X[i,]=c(g1(temp[1]),g2(temp[2]),g3(temp[3]),g4(temp[4]),g5(temp[5]))
   cat(X[i,],"\n")
 }

 X 


 t=seq(0.01,10,length=2000)
 plot(t,t+log(t),type="l",col=2)
 g1<-function(x) (x+exp(-x))/2
 g2<-function(x) exp(-x)
 g3<-function(x) -log(x)
 
 Loop=20
 X=matrix(0.1,Loop,3)
 for(i in 2:Loop){
   temp=X[i-1,]
   X[i,]=c(g1(temp[1]),g2(temp[2]),g3(temp[3]))
   cat(X[i,],"\n")
 }

 X 

G1<-function(x) (1+1/x-log(x))/(1+x)+x
G2<-function(x) 2*(1+1/x-log(x))/(1+x)+x
G3<-function(x) 3*(1+1/x-log(x))/(1+x)+x
G4<-function(x) 4*(1+1/x-log(x))/(1+x)+x
G5<-function(x) 5*(1+1/x-log(x))/(1+x)+x

 Loop=200
 X=matrix(0.21,Loop,5)
 for(i in 2:Loop){
   temp=X[i-1,]
   X[i,]=c(G1(temp[1]),G2(temp[2]),G3(temp[3]),G4(temp[4]),G5(temp[5]))
   cat(X[i,],"\n")
 }

 X 


  logistic<-function(data1,delta=0.0001,Loop=100){
    X=data1$X
    Y=data1$Y
    tol = .Machine$double.eps^0.5
    loop=0
    p=ncol(X)
    n=nrow(X)
    beta_old=rep(1,p)
    beta_new=rep(0,p)
    while((sum((beta_new-beta_old)^2)>1e-13)&loop<Loop){
      loop=loop+1
      beta_old=beta_new      
      temp=1/(1+exp(X%*%beta_old))
      Wt=diag(as.vector(temp))
      Sn=-t(X)%*%Y+t(X)%*%temp #apply(t(X)%*%Wt,1,sum)   
      Temp=temp/(1+exp(-X%*%beta_old))
      WT=diag(as.vector(Temp))
      Dn=-t(X)%*%WT%*%X      
      beta_new=beta_old-solve(Dn+delta*diag(p)/n)%*%Sn      
    } 
    list(beta.hat=beta_new,hessian=Dn,loop=loop)
  }  


  #鸢尾花
  n=100
  u=runif(n)
  beta=c(1,-1,1,-1,0)
  p=length(beta)
  X=matrix(rnorm(n*p),n,p)
  temp=exp(X%*%beta)+1
  p0=1/temp
  Y=1*(u<p0)
  X=cbind(1,X)
  data1=list(X=X,Y=Y)
  fit=logistic(data1,delta=1)
  fit$beta.hat
  solve(-fit$hessian)*n







#Iris Data Set（鸢尾属植物数据集）是历史最悠久的数据集，
#它首次出现在著名的英国统计学家和生物学家Ronald Fisher 
#1936年的论文《The use of multiple measurements in taxonomic problems》中，
#被用来介绍线性判别式分析。在这个数据集中，包括了三类不同的鸢尾属植物：
#Iris Setosa，Iris Versicolour，Iris Virginica。
#每类收集了50个观测值，因此这个数据集一共包含了150个观测值。
#
###该数据集测量了所有150个样本的4个特征，分别是：
#1. sepal length（花萼长度）
#2. sepal width（花萼宽度）
#3. petal length（花瓣长度）
#4. petal width（花瓣宽度）
#以上四个特征的单位都是厘米（cm）
 source("logistic.r")
 mydata=iris
 mydata1=mydata[!mydata[,5]=="setosa",]
 
 #mydata1=subset(iris,Species=="virginica"|Species=="versicolor")
 fit1<-glm(Species~.,data=mydata1,family="binomial")
 summary(fit1)
 X=mydata1[,1:4]
 X=as.matrix(X)
 Y=1-1*(mydata1[,5]=="virginica")
 mydata2=list(X=cbind(1,X),Y=Y)
 fit2=logistic(mydata2)
 beta.hat=fit2$beta.hat
 beta.sd=sqrt(diag(solve(-fit2$hessian)))
 cbind(beta.hat,beta.sd,fit1$coef, sqrt(diag(vcov(fit1))))
 
 






######################牛顿迭代
#等高线（contour map)是可视化二维空间标量场的基本方法，可以将
#三维数据使用二维的方法可视化，同时用颜色视觉特征表示第三维数据，如地图上的等高线、天气预报中的
#等压线和等温线等等。假设f(x,y)在点(x,y)处的数值，等高线是在二维数据场中满足f(x,y)=c的空间点集合
#按一定的顺序连接而成的线。数值为c的值线可以将二维空间坐标场分为两部分：如果f(x,y)<c,则该点在等值线内；
#如果f(x,y)>c,则该点在等值线外
 library(reshape2)
 library(ggplot2)
 library(directlabels) 
 library(RColorBrewer)

 


 myf<-function(x,y) (1-x)^2+100*(y-x^2)^2+0.3*(0.2-2*y)^2+100*(x-y^2)^2-0.5*(x^2+5*y^2)
 myg<-function(x,y) c(2*(x-1)+400*x*(x^2-y)+200*(x-y^2)-x,200*(y-x^2)+1.2*(2*y-0.2)+400*(y^2-x)*y-5*y)
 myG<-function(x,y) matrix( c(201+400*(x^2-y)+800*x^2,-400*(x+y),-400*(x+y),197.4+1200*y^2-400*x),2,2)


 n=100
 x <- seq(-1, 2, length.out = n)
 y <- seq(-1, 2, length.out = n)
 z<-outer(x,y,myf)
 persp(x,y,z,theta=0,phi=15,expand=0.7,col="lightblue",xlab="X",ylab="Y",zlab="Z")

  m=6
  RES1=matrix(1.5,m,2)
  RES1[1,]=c(0.5,-1)
  for(i in 2:m){
    RES1[i,]=RES1[i-1,]-solve(myG(RES1[i-1,1],RES1[i-1,2]))%*%myg(RES1[i-1,1],RES1[i-1,2])
  }
 
  round(RES1,4)


  m2=7
  RES2=matrix(0,m2,2)
  RES2[1,]=c(1.5,-1)
  for(i in 2:m2){
    RES2[i,]=RES2[i-1,]-solve(myG(RES2[i-1,1],RES2[i-1,2]))%*%myg(RES2[i-1,1],RES2[i-1,2])
  }
   round(RES2,4)


   m3=8
  RES3=matrix(0,m3,2)
  RES3[1,]=c(1.5,1.5)
  for(i in 2:m3){
    RES3[i,]=RES3[i-1,]-solve(myG(RES3[i-1,1],RES3[i-1,2]))%*%myg(RES3[i-1,1],RES3[i-1,2])
  }
   round(RES3,4) 
  contour(x, y, log10(z+2.1),xlab = "X", ylab = "Y",col=5, nlevels = 17)

  points(RES1[1,1], RES1[1,2], cex = 2, pch = 20,col=3)  
  points(RES1[m,1], RES1[m,2], cex = 1, pch = 20,col=4)
  lines(RES1, pch=3, type="o",col=2)

  points(RES2[1,1], RES2[1,2], cex = 2, pch = 20,col=3)
  points(RES2[m2,1], RES2[m2,2], cex = 1, pch = 20,col=4)
  lines(RES2, pch=3, type="o",col=2)

  points(RES3[1,1], RES3[1,2], cex = 2, pch = 20,col=3)
  points(RES3[m3,1], RES3[m3,2], cex = 1, pch = 20,col=4)
  lines(RES3, pch=3, type="o",col=2)
  c(myf(RES1[m,1],RES1[m,2]),myf(RES2[m2,1],RES2[m2,2]),myf(RES3[m3,1],RES3[m3,2]))


##############################



 rm(list = ls())
 myf<-function(x,y) (1-x)^2+100*(y-x^2)^2+0.3*(0.2-2*y)^2+100*(x-y^2)^2-0.5*(x^2+5*y^2)
 myg<-function(x,y) c(2*(x-1)+400*x*(x^2-y)+200*(x-y^2)-x,200*(y-x^2)+1.2*(2*y-0.2)+400*(y^2-x)*y-5*y)
 myG<-function(x,y) matrix( c(201+400*(x^2-y)+800*x^2,-400*(x+y),-400*(x+y),197.4+1200*y^2-400*x),2,2)
 x3=c(1.8,1.28) #初值
 x2=c(1.5,-1)  #初值
 x1=c(0.5,-1)  #初值
 n=100
 x <- seq(-1.2, 2.2, length.out = n)
 y <- seq(-1.2, 2.2, length.out = n)
 z<-outer(x,y,myf)
 #persp(x,y,z,theta=0,phi=15,expand=0.7,col="lightblue",xlab="X",ylab="Y",zlab="Z")
 m=6
 RES1=matrix(1.5,m,2)
 RES1[1,]=x1
 for(i in 2:m){
   RES1[i,]=RES1[i-1,]-solve(myG(RES1[i-1,1],RES1[i-1,2]))%*%myg(RES1[i-1,1],RES1[i-1,2])
 } 
 round(RES1,4) 
 m2=7
 RES2=matrix(0,m2,2)

 RES2[1,]=x2
 for(i in 2:m2){
   RES2[i,]=RES2[i-1,]-solve(myG(RES2[i-1,1],RES2[i-1,2]))%*%myg(RES2[i-1,1],RES2[i-1,2])
 }
 round(RES2,4)

 m3=6
 RES3=matrix(0,m3,2)
 
 RES3[1,]=x3
 for(i in 2:m3){
   RES3[i,]=RES3[i-1,]-solve(myG(RES3[i-1,1],RES3[i-1,2]))%*%myg(RES3[i-1,1],RES3[i-1,2])
 }
  round(RES3,4) 
 contour(x, y, log10(z+2.1),xlab = "X", ylab = "Y",col=5, nlevels = 17)




  image(x, y, log10(z+2.1), main="图8",col = terrain.colors(200), xlab = "X", ylab = "Y")
  contour(x, y,log10(z+2.1),col=5, nlevels = 17, add=T)
  points(RES1[1,1], RES1[1,2], cex = 2, pch = 20,col=3)  
  points(RES1[m,1], RES1[m,2], cex = 1, pch = 20,col=4)
  #lines(RES1, pch=3, type="o",col=2)

  s <- seq(dim(RES1)[1]-1)
  arrows(RES1[, 1][s], RES1[, 2][s], RES1[, 1][s+1], RES1[, 2][s+1], col=2, length=0.1, lwd=2.5)

  points(RES2[1,1], RES2[1,2], cex = 2, pch = 20,col=3)
  points(RES2[m2,1], RES2[m2,2], cex = 1, pch = 20,col=4)
  #lines(RES2, pch=3, type="o",col=2)

  s <- seq(dim(RES2)[1]-1)
  arrows(RES2[, 1][s], RES2[, 2][s], RES2[, 1][s+1], RES2[, 2][s+1], col=2, length=0.1, lwd=2.5)

  points(RES3[1,1], RES3[1,2], cex = 2, pch = 20,col=3)
  points(RES3[m3,1], RES3[m3,2], cex = 1, pch = 20,col=4)
  #lines(RES3, pch=3, type="o",col=2)

  s <- seq(dim(RES3)[1]-1)
  arrows(RES3[, 1][s], RES3[, 2][s], RES3[, 1][s+1], RES3[, 2][s+1], col=2, length=0.1, lwd=2.5)


  c(myf(RES1[m,1],RES1[m,2]),myf(RES2[m2,1],RES2[m2,2]),myf(RES3[m3,1],RES3[m3,2]))


##############################
##################################################33

  myf<-function(x,y) log10(3*x^2+3*y^2-x^2*y+1)

   myf1<-function(x,y)  (3*x^2+3*y^2-x^2*y+1)

  myg<-function(x,y) c(6*x-2*x*y,6*y-x^2)
  myG<-function(x,y) matrix(c(6-2*y,-2*x,-2*x,6),2,2)

  n=100
  x <- seq(-6, 6, length.out = n)
  y <- seq(-6, 6, length.out = n)
  z<-outer(x,y,myf)
  z1<-outer(x,y,myf1)
  persp(x,y,z1,theta=0,phi=15,expand=0.7,col="lightblue",xlab="X",ylab="Y",zlab="Z")

  par(mar = c(4,4,0.5,0.5))
  contour(x, y, z,xlab = "x", ylab = "y", nlevels = 10)
  m=7
  RES1=matrix(1.5,m,2)
  for(i in 2:m){
    RES1[i,]=RES1[i-1,]-solve(myG(RES1[i-1,1],RES1[i-1,2]))%*%myg(RES1[i-1,1],RES1[i-1,2])
  }
 
  round(RES1,4)


  m=6
  RES2=matrix(0,m,2)
  RES2[1,]=c(2,4)
  for(i in 2:m){
    RES2[i,]=RES2[i-1,]-solve(myG(RES2[i-1,1],RES2[i-1,2]))%*%myg(RES2[i-1,1],RES2[i-1,2])
  }
   round(RES2,4)

  par(mfrow = c(1, 2))
  par(mar = c(4,4,0.5,0.5))
  contour(x, y, z,xlab = "x", ylab = "y", nlevels = 10)
  lines(RES1, pch=3, type="o",col=2)  
  contour(x, y, z,xlab = "x", ylab = "y", nlevels = 10)
  lines(RES2, pch=4, type="o",col=4)
  
  #取初值为(0,3)牛顿方法失效，因为Hessian 矩阵不可逆
   myG(0,3)





########有问题，学生写的程序，没有时间修改

library(ggplot2)
library(ggalt)
library(ggforce)
k=6.5
p=1
r=1
tau=0.5
sigma=1
x=seq(-5,5,len=500)

Y1=function(x) k*dnorm(x)
y1<-Y1(x)
Y2<-function(x) (tau*sigma^2)^(-r*p-p/2)*exp(-x^2/2/tau/sigma^2)*x^(2*r)
y2<-Y2(x)
mu=seq(0,0,length.out = 500)
mu0=seq(0,0,length.out = 500)
for (i in 1:500) {
  mu[i]<-sample(runif(x,0,y1[i]),size=1) 

  if(mu[i]<y2[i]) mu0[i]=mu[i]
  if(mu0[i]==0) mu0[i]=NA 
}

data<-as.data.frame(cbind(x,y1,y2))


dt<-as.data.frame((cbind(data,mu,mu0)))
ggplot(dt,aes(x,))+
  geom_line(aes(,y1),color='blue',size=1.5)+
 
  geom_line(aes(,y2),color='red',size=1.5)+
  
  geom_ribbon(aes(x=x, ymin=y2, ymax = y1),
              alpha = 0.3)+
geom_point(aes(,mu),color='green',size=1)+
  geom_point(aes(,mu0),color='black',size=1)+
geom_segment(aes(x=1.06060606,y=Y1(1.06060606),xend=1.06060606,yend=Y2(1.06060606)),colour='yellow'
             ,lwd=1.3)+
  geom_segment(aes(x=1.06060606,y=Y2(1.06060606),xend=1.06060606,yend=0),colour='purple'
               ,lwd=1.3)+
  geom_hline(yintercept = 0)
###########################