#A Kernel Method for the Two-Sample-Problem##d.K都是单个写的
rm(list=ls())
set.seed(123)
library(MASS)
library(Rcpp)
library(PMMD)
t1 <- Sys.time()
n<-100
m<-100
d<-20#c(20,40,80,100)

B=100
Sim=100
alpha=0.05
ECp <- matrix(0,Sim,B)
KSECp <- matrix(0,Sim,B)
#p_value<-rep(0,Sim)

for (sim in 1:Sim) {
  X <- matrix(rnorm(n*d, 0, 1), n, d)
  Y <- matrix(rnorm(m*d, 0.05, 1), m, d)

  Tn<-meammd(X=X, Y=Y, pval=TRUE, type="proj", numproj=10,numperm = 0,seednum = 123)$stat
  KSTn <- sqrt(maxmmd(X=X, Y=Y, pval=TRUE, type="proj", numproj=10,numperm = 0,seednum = 123))
  for (b in 1:B) {
    cat(sim,b,"\r")
    Z<-rbind(X,Y)
    set<-sample(1:(n+m),n)
    a<-as.vector(set)

    X1<-as.matrix(Z[a,])
    Y1<-as.matrix(Z[-a,])

    Tn_b<-meammd(X=X1, Y=Y1, pval=TRUE, type="proj", numproj=10,numperm =0,seednum = 123)$stat
    KSTn_b<-sqrt(maxmmd(X=X1, Y=Y1, pval=TRUE, type="proj", numproj=10,numperm =0,seednum = 123))
    ECp[sim,b]<-1*(Tn_b>=Tn)
    KSECp[sim, b] <- 1*(KSTn_b>=KSTn)
  }
}
ECP<-sum(1*((rowMeans(ECp))<0.05))/Sim
KSECP <- sum(1*((rowMeans(KSECp))<0.05))/Sim
t2 <- Sys.time()
print((t2-t1))
print(ECP)
print(KSECP)

