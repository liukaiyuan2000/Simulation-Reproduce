rm(list = ls())
library(kernlab)
library(energy)
source("D:/Desktop/JASAIndependenceTest/8014586/source.R")
source("D:/Desktop/JASAIndependenceTest/code/DGP.R")
Rcpp::sourceCpp("D:/Desktop/JASAIndependenceTest/code/function.cpp")
B = 50
nsim = 1000
RES = NULL
for(n in c(25, 50, 100)){
  for(model.choose in 1:4){
    res.slice  = rep(0, nsim)
    res.kernel = rep(0, nsim)
    res.dcov   = rep(0, nsim)
    res.disco  = rep(0, nsim)
    for(i in 1:nsim){
      cat('Example 1: ', 'n = ', n, '; model = ', model.choose, '; ', i, '\r', sep = "")
      dat = DGP(example = 1, model = model.choose, n = n, p = 5, q = 1)
      x = dat$x
      y = dat$y
      y.d = dat$y.d
      res.slice[i] = ptest_slm_rcpp(x, y.d, B = B)$pvalue
      res.kernel[i] = ptest_kem_rcpp(x, y, B = B)$pvalue
      res.dcov[i] = dcov.test(x, y, R = B)$p
      res.disco[i] = disco(x, y.d, R = B)$p
    }
    cat('\n  Done. \n')
    temp = cbind(res.slice, res.kernel, res.dcov, res.disco)
    RES = rbind(RES, sprintf('%.4f', colMeans(temp <= 0.1)))
  }
}
RES = cbind(
  c('(a)', NA, NA, '(b)', NA, NA,
    '(c)', NA, NA, '(d)', NA, NA),
  rep(c(25, 50, 100), 4), RES
)
colnames(RES) = c('Model', 'n', 'Slice', 'Kernel', 'DCOV', 'DISCO')
as.data.frame(RES)

save.image("D:/Desktop/JASAIndependenceTest/results/Example 1.RData")

library(stargazer)
stargazer(RES, summary = F, digits = 4, table.placement = "H")

