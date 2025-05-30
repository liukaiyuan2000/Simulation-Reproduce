setwd('D:/Desktop/mean test')
Rcpp::sourceCpp('CLX.cpp')
Rcpp::sourceCpp('CQ.cpp')
Rcpp::sourceCpp('SDu.cpp')
Rcpp::sourceCpp('aSPUperm.cpp')
library(highmean)
library(MASS)
set.seed(11)

nsim = 2000
n1 <- n2 <- 50
p <- 200
B <- 300
true.cov <- 0.6^(abs(outer(1:p, 1:p, "-"))) # AR1 covariance

DGP <- function(r){
  mu.cons <- (2*r*(1/n1+1/n2)*log(p))^0.5
  mu1 <- rep(0, p)
  if(r == 0){
    mu2 <- mu1
  }else{
    mu2[sample(1:p, 117, replace = F)] = mu.cons
  }
  sam1 <- mvrnorm(n = n1, mu = mu1, Sigma = true.cov)
  sam2 <- mvrnorm(n = n2, mu = mu2, Sigma = true.cov)
  return(
    list(sam1, sam2)
  )
}

t1 = Sys.time()
RES = NULL
for(r in seq(0, 0.08, 0.02)){
  res = matrix(0, nsim, 13)
  for(i in 1:nsim){
    cat('r = ', r, '; ', i, '\r', sep = '')
    dat <- DGP(r)
    sam1 <- dat[[1]]
    sam2 <- dat[[2]]
    res.SPU = aSPUperm_rcpp(sam1, sam2, n_perm = B)
    res.CLZ = epval_Chen2014(sam1, sam2, n.iter = B)$pval
    res.CLX = CLX(sam1, sam2)
    res.BS  = epval_Bai1996(sam1, sam2, n.iter = B)$pval
    res.CQ  = CQ(sam1, sam2)
    res.SD  = SDu(sam1, sam2)
    res[i, ] = c(
      res.SPU, res.CLZ, res.CLX, res.BS, res.CQ, res.SD
    )
  }
  cat('\n   Done! \n')
  RES = cbind(RES, 100*colMeans(res <= 0.05))
}
t2 = Sys.time()
print(difftime(t2, t1, units = 'mins'))
rnames = c(
  'SPU1', 'SPU2', 'SPU3', 'SPU4', 'SPU5', 'SPU6', 
  'SPU_Inf', 'aSPU', 'CLZ', 'CLX', 'BS', 'SQ', 'SD'
)
cnames = c(
  'r = 0', 'r = 0.02', 'r = 0.04', 
  'r = 0.06', 'r = 0.08'
)
RES.dat = as.data.frame(
  RES, row.names = rnames
)
colnames(RES.dat) = cnames
round(RES.dat)

