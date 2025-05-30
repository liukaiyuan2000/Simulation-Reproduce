#Beijing University of Technology
#6 May 2022
#Statistical methods for linear quantile models with a censored covariate
rm(list =ls())
library(MASS)
library(survival)
library(quantreg)
source("BTF.r")
tau <- 0.25
beta0 <- c(-2, 2, 1)
p = length(beta0);
rho = 0.5
Sigmma = matrix(0, p, p)
Sigmma <- (rho) ^ abs(row(Sigmma) - col(Sigmma)) 
N <- 200 * (1:2) #样本量
Sim <- 1000 #模拟次数
RES <- NULL
#执行循环----
for(n in N){
  BETAN <- matrix(0, Sim, p)
  BETAM <- matrix(0, Sim, p)
  for(sim in 1:Sim){
    cat(n, sim, "\r")
    e <- rnorm(n)
    e <- e - quantile(e, tau)
    X <- mvrnorm(n, rep(0, p), Sigmma) #X ~ N(0, Sigma)
    Y <- X %*% beta0 + e # ture model
    C <- rexp(n, 1) - 1 #删失水平
    Xc <- pmin(X[, 1], C) #观测值, X第一个分量发生删失
    delta <- 1 * (X[, 1] < C) #若为非删失数据则delta为1
    mean(delta) #删失率
    W <- BTF(Xc, delta)
    W <- drop(W)
    XX <- cbind(Xc, X[, -1])
    beta.hat <- rq(Y ~ XX + 0, tau = tau, weights = W)$coef
    beta.n <- rq(Y ~ XX + 0, tau = tau)$coef
    BETAN[sim,] <- beta.n - beta0
    BETAM[sim,] <- beta.hat - beta0
  }
  resm <- rbind(colMeans(BETAN), colMeans(BETAM))
  ress <- rbind(apply(BETAN, 2, sd), apply(BETAM, 2, sd))
  RES <- rbind(RES, cbind(resm, ress))
}

#输出结果----
#前两行为n=200,后两行n=400
print(round(RES * 100, 2))











