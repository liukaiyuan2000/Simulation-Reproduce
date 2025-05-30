x=seq(0,5,length=1000)
g <- function(x) log(x)/(1+x)
g1 <- function(x) (1+1/x-log(x))/(1+x)/(1+x)
root1 <- uniroot(g1,interval = c(2,4))$root
main1 <- expression(x %~~% 3.59112)
plot(x,g(x),ylim=c(0,0.3),type="l",lwd=2, main="ͼ1")
abline(v=root1,col=2,lwd=2,lty=2)
points(root1,g(root1),pch=19)
text(3.2,0.1,main1)

t_star <- expression(x^"*")
t0 <- expression(x^(0))
t1 <- expression(x^(1))
t2 <- expression(x^(2))
b0 <- expression(b[0])
b1 <- expression(b[1])
b2 <- expression(b[2])
a0 <- expression(a[0])
a1 <- expression(a[1])
a2 <- expression(a[2])

x0 <- c(1,5)
y0 <- c(-0.1,-0.1)
x1 <- c(3,5)
y1 <- c(-0.2,-0.2)
x2 <- c(3,4)
y2 <- c(-0.3,-0.3)
plot(x,g1(x),ylim=c(-0.4,0.5),xlim=c(0.5,5.5),family="serif",ylab="g'(x)",type="l",lwd=2, main="ͼ2")
abline(h=0,col=2,lwd=1,lty=2)
lines(x0,y0,xlim=c(1,5),lwd=2)
points(1,-0.1,pch=19)
points(3,-0.1,pch=19)
points(5,-0.1,pch=19)
lines(x1,y1,xlim=c(3,5),lwd=2)
points(3,-0.2,pch=19)
points(4,-0.2,pch=19)
points(5,-0.2,pch=19)
lines(x2,y2,xlim=c(4,5),lwd=2)
points(3,-0.3,pch=19)
points(3.5,-0.3,pch=19)
points(4,-0.3,pch=19)
points(root1,g1(root1),pch=19)
text(0.7,-0.1,a0)
text(2.7,-0.2,a1)
text(2.7,-0.3,a2)
text(3,-0.03,t0)
text(4,-0.13,t1)
text(3.5,-0.23,t2)
text(5.3,-0.1,b0)
text(5.3,-0.2,b1)
text(4.3,-0.3,b2)
text(root1,g1(root1)+0.07,t_star)

aa <- 3-g1(3)/g2(3)
x0 <- c(3-g1(3)/g2(3),3-g1(3)/g2(3))
y0 <- c(0,g1(3-g1(3)/g2(3)))
g1 <- function(x) (1+1/x-log(x))/(1+x)/(1+x)
g2 <- function(x) ((-1/x-1/x^2)*(1+x)^2-2*(1+x)*(1/x+1-log(x)))/(1+x)^4
plot(x,g1(x),xlim=c(2.6,3.9),ylim=c(-0.01,0.03),family="serif",ylab="g'(x)",type="l",lwd=2, main="ͼ3")
abline(h=0,col=2,lwd=1,lty=2)
points(3,g1(3),pch=19)
lines(x,g1(3)+g2(3)*(x-3),col=3,lty=5,lwd=1)
lines(x0,y0,col=3,lty=2,lwd=1)
points(3-g1(3)/g2(3),g1(3-g1(3)/g2(3)),pch=19)
points(aa-g1(aa)/g2(aa),g1(aa-g1(aa)/g2(aa)),pch=19)
points(root1,g1(root1),pch=19)
text(3.08,g1(3)+0.002,t0)
text(3.44,0.007,t1)
text(aa-g1(aa)/g2(aa),0.0038,t2)
text(root1,-0.003,t_star)



f<-function(x) 0.25*(x-1)^4+exp(x)
f1<-function(x) (x-1)^3+exp(x)
f2<-function(x) 3*(x-1)^2+exp(x)
q<-function(x,x0) f(x0)+f1(x0)*(x-x0)+0.5*f2(x0)*(x-x0)^2
t=seq(-1,3.5,length=2000)
x0=2.5
t1=seq(1,3,length=200)  
plot(t,f(t),type='l',ylim=c(-5,45),xlab="x",ylab="f(x)",xlim=c(-1,4),bty="L",lwd=2)
par(new=TRUE)
plot(t1, q(t1,x0),type='l',ylim=c(-5,45),xlab="x",ylab="f(x)",col=2,xlim=c(-1,4),lwd=2)


x1=x0-f1(x0)/f2(x0) 
t2=seq(x1-1.2,x1+0.5,length=200)
par(new=TRUE)
plot(t2, q(t2,x1),lwd=2,type='l',ylim=c(-5,45),xlab="x",ylab="f(x)",col=3,xlim=c(-1,4))

x2=x1-f1(x1)/f2(x1)
par(new=TRUE)
t3=seq(x2-1,x2+1,length=200) 
plot(t3, q(t3,x3),lwd=2,type='l',ylim=c(-5,45),xlab="x",ylab="f(x)",col=6,xlim=c(-1,4))


x3=x2-f1(x2)/f2(x2)
lines(rep(x0,200),seq(-0.10,f(x0),length=200),lty=2,col=4,lwd=1)
lines(rep(x1,200),seq(-0.10,f(x1),length=200),lty=2,col=4,lwd=1)
lines(rep(x2,200),seq(-0.10,f(x2),length=200),lty=2,col=4,lwd=1)
lines(rep(x3,200),seq(-0.10,f(x3),length=200),lty=2,col=4,lwd=1)
abline(h=0,col=2,lwd=1,lty=2)

text(x0,-2,expression(x[0]==2.5))
text(x1,-2,expression(x[1]))
text(x2,-2,expression(x[2]))
text(x3,-2,expression(x[3]))
###############################
###############################
###############################


x0 <- c(1.5-g1(1.5)/(g1(2.5)-g1(1.5)),1.5-g1(1.5)/(g1(2.5)-g1(1.5)))
y0 <- c(0,g1(1.5-g1(1.5)/(g1(2.5)-g1(1.5))))
t_star <- expression(x^"*")
t0 <- expression(x^(0))
t1 <- expression(x^(1))
t2 <- expression(x^(2))
t3 <- expression(x^(3))
plot(x,g1(x),ylim=c(-0.05,0.23),xlim=c(1.3,4),family="serif",ylab="g'(x)",type="l",lwd=2, main="ͼ5")
points(root1,g1(root1),pch=19)
points(1.5,g1(1.5),pch=19)
points(2.5,g1(2.5),pch=19)
abline(h=0,col=2,lwd=2,lty=2)
lines(x,g1(1.5)+(g1(2.5)-g1(1.5))*(x-1.5),lty=2,lwd=2)
points(1.5-g1(1.5)/(g1(2.5)-g1(1.5)),g1(1.5-g1(1.5)/(g1(2.5)-g1(1.5))),pch=19)
lines(x0,y0,lty=2,col=3,lwd=2)
points(ab,g1(ab),pch=19)
text(1.65,0.22,t0)
text(2.5,0.07,t1)
text(2.75,0.05,t2)
text(3.2,0.03,t3)
text(root1,0.02,t_star)
#################################







f<-function(t) log(t)/(1+t)
f1<-function(t) (1+1/t-log(t))/(1+t) 
x=seq(1.5,5,length=1000)

X=matrix(0,10,3)
X[1,2:3]=c(1.5,2)
for(i in 2:10){   
  temp=X[i-1,3]-f1(X[i-1,3])*(X[i-1,3]-X[i-1,2])/(f1(X[i-1,3])-f1(X[i-1,2]))
  X[i,]=c(X[i-1,2],X[i-1,3],temp)
} 

x1 <- c(1.5, 5)
y1 <- c(0,0)
plot(x,f1(x),type="l",col=1,lwd=1.5,ylab=expression(f(x)),main="ͼ6")
abline(h=0,col=2,lwd=1,lty=2)
for(i in 2:6){
  points(X[i, 1], f1(X[i, 1]), pch=19)
}

for(i in 2:5){
  xt <- X[i, 3]+0.1
  yt <- f1(X[i, 1])+(f1(X[i+1, 1])-f1(X[i, 1]))/(X[i+1, 1]-X[i, 1])*(xt-X[i, 1])
  xtt <- X[i, 1]-0.2
  ytt <- f1(X[i, 1])+(f1(X[i+1, 1])-f1(X[i, 1]))/(X[i+1, 1]-X[i, 1])*(xtt-X[i, 1])
  x1 <- c(xtt, xt)
  y1 <- c(ytt, yt)
  lines(x1, y1, lty=i+1, lwd=3,col=i+1)
}

for(i in 2:6){
  x1 <- c(X[i, 1], X[i, 1])
  y1 <- c(f1(X[i, 1]), 0)
  lines(x1, y1, lty=2, lwd=3)
}

text(X[2, 1], -0.02, expression(x^{(0)}))
text(X[2, 2], -0.02, expression(x^{(1)}))
text(X[2, 3], -0.02, expression(x^{(2)}))
text(X[3, 3], -0.02, expression(x^{(3)}))
text(X[4, 3], -0.02, expression(x^{(4)}))



aa <- 1.75+4*g1(1.75)
ab <- aa+4*g1(aa)
ac <- ab+4*g1(ab)
x0 <- c(aa,aa)
y0 <- c(0,g1(aa))
x1 <- c(ab,ab)
y1 <- c(0,g1(ab))
plot(x,g1(x),ylim=c(-0.05,0.23),xlim=c(1.3,4),family="serif",ylab="g'(x)",type="l",lwd=2, main="ͼ7")
abline(h=0,col=2,lwd=1,lty=2)
points(1.75,g1(1.75),pch=19)
points(aa,g1(aa),pch=19)
points(ab,g1(ab),pch=19)
points(ac,g1(ac),pch=19)
points(root1,g1(root1),pch=19)
lines(x,g1(1.75)-0.25*(x-1.75),lty=2,col=3)
lines(x0,y0,lty=3,lwd=2)
lines(x,g1(aa)-0.25*(x-aa),lty=2,col=4)
lines(x1,y1,lty=3,lwd=2)
text(1.7,0.12,t0)
text(2.35,0.08,t1)
text(2.6,0.06,t2)
text(2.75,0.05,t3)
text(root1,0.02,t_star)


