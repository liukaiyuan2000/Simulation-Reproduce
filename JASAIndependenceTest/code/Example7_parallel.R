rm(list = ls())
library(ECCFIC2024)
library(doSNOW)
library(tcltk)
no_cores <- 12
progress <- function(n){
  pb <- txtProgressBar(max = nsim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)
B = 50
nsim = 1000
RES = NULL
cl <- makeSOCKcluster(no_cores)
registerDoSNOW(cl)
t1 <- Sys.time()
for(n in c(10, 20, 50)){
  for(model.choose in 1:6){
    cat('Example 7: ', 'n = ', n, '; model = ', model.choose, '\n', sep = "")
    temp = foreach(
      sim = 1:nsim, .options.snow = opts, 
      .packages = c('ECCFIC2024', 'energy'), .combine = 'rbind'
    ) %dopar% {
      dat = DGP(example = 7, model = model.choose, n = n)
      x = dat$x
      y = dat$y
      y.d = dat$y.d
      return(c(
        ptest_slm_rcpp(x, y.d, B = B)$pvalue, 
        ptest_kem_rcpp(x, y, B = B)$pvalue, 
        dcov.test(x, y, R = B)$p, 
        disco(x, y.d, R = B)$p
      ))
    }
    cat('\n  Done. \n')
    RES = rbind(RES, sprintf('%.4f', colMeans(temp <= 0.1)))
  }
}
t2 <- Sys.time()
stopCluster(cl)
print((t2-t1))
RES = cbind(
  c('A(1)', NA, NA, 'A(2)', NA, NA, 
    'A(3)', NA, NA, 'B(1)', NA, NA, 
    'B(2)', NA, NA, 'B(3)', NA, NA), 
  rep(c(10, 20, 50), 6), RES
)
colnames(RES) = c('Model', 'n', 'Slice', 'Kernel', 'DCOV', 'DISCO')
as.data.frame(RES)

save.image("D:/Desktop/JASAIndependenceTest/results/Example 7.RData")

library(stargazer)
stargazer(RES, summary = F, digits = 4, table.placement = "H")

