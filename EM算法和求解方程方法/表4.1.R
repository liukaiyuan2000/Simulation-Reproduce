x = c(85, 196, 341)
n = rep(0,6)
p = rep(1/3,3)
itr = 40
## EXPECTATION AND MAXIMIZATION FUNCTIONS
allele.e = function(x,p){
  n.cc = (x[1]*(p[1]^2))/((p[1]^2)+2*p[1]*p[2]+2*p[1]*p[3])
  n.ci = (2*x[1]*p[1]*p[2])/((p[1]^2)+2*p[1]*p[2]+2*p[1]*p[3])
  n.ct = (2*x[1]*p[1]*p[3])/((p[1]^2)+2*p[1]*p[2]+2*p[1]*p[3])
  n.ii = (x[2]*(p[2]^2))/((p[2]^2)+2*p[2]*p[3])
  n.it = (2*x[2]*p[2]*p[3])/((p[2]^2)+2*p[2]*p[3])
  n = c(n.cc,n.ci,n.ct,n.ii,n.it,x[3])
  return(n)
}

allele.m = function(x,n){
  p.c = (2*n[1]+n[2]+n[3])/(2*sum(x))
  p.i = (2*n[4]+n[5]+n[2])/(2*sum(x))
  p.t = (2*n[6]+n[3]+n[5])/(2*sum(x))
  p = c(p.c,p.i,p.t)
  return(p)
}


for(i in 1:itr){
  n = allele.e(x,p)
  p = allele.m(x,n)
}
round(p,5)


p2 <- matrix(0,nrow=9,ncol=3)
p2[1,] <- rep(1/3,3)


for(i in 1:8){
  n = allele.e(x,p2[i,])
  p2[(i+1),] = allele.m(x,n)
}
p_t <- p2[,1:2]
round(p_t,6)


R <- 0
for(i in 1:8){
  R[i]=sqrt((p_t[i+1,1]-p_t[i,1])^2+(p_t[i+1,2]-p_t[i,2])^2)/sqrt(p_t[i,1]^2+p_t[i,2]^2)
}
round(R,7)


D_c <- 0
for(i in 1:8){
  D_c[i]=(p_t[i+1,1]-p[1])/(p_t[i,1]-p[1])
}
round(D_c,4)


D_i <- 0
for(i in 1:8){
  D_i[i]=(p_t[i+1,2]-p[2])/(p_t[i,2]-p[2])
}
round(D_i,3)


b4_1 <- cbind(round(p_t,6),round(c(NA,R),7),round(c(NA,D_c),4),round(c(NA,D_i),3))

xnames <- seq(0,8,by=1)
ynames <- c(expression(p[C]^(t)),expression(p[I]^(t))
            ,expression(R^(t)),expression(D[C]^(t)),
            expression(D[I]^(t)))
b4_1 <- matrix(b4_1,nrow=9,dimnames=list(xnames,ynames))
