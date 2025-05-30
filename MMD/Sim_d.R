# A Kernel Method for the Two-Sample-Problem
rm(list = ls())
library(PMMD)
library(doSNOW)
library(tcltk)
library(stargazer)
setwd('D:/Desktop/MMD')
############# SETTING AND FUNCTION ##################
d.choose <- seq(10, 400, by = 20)
n = 50
m = 50
B = 100
Sim = 500
alpha = 0.05
K = 10

no_cores <- 12
progress <- function(n){
  pb <- txtProgressBar(max = Sim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)
#######################################################

##################### Simulation ######################
t1 <- Sys.time()
cl <- makeSOCKcluster(no_cores)
registerDoSNOW(cl)
RES.cvm = c()
for(d in d.choose){
  cat("n = ", n, "; d = ", d, "\n", sep = "")
  res <- foreach(
    sim = 1:Sim, .options.snow = opts, 
    .packages = c('PMMD', 'MASS'), .combine = 'rbind'
  )  %dopar% {
    Sigma = 0.2^abs(outer(1:d, 1:d, '-'))
    X <- matrix(rnorm(n*d), n, d)
    Y <- mvrnorm(n, rep(0, d), Sigma)
    
    Tn <- meammd(
      X = X, Y = Y, pval = TRUE, type = "proj", numproj = K,
      numperm = 0
    )$stat
    Tn_b = rep(0, B)
    for (b in 1:B) {
      Z <- rbind(X, Y)
      set <- sample(1:(n+m), n)
      
      X1 <- Z[set, ]
      Y1 <- Z[-set,]
      
      Tn_b[b] <- meammd(
        X = X1, Y = Y1, pval = TRUE, type = "proj", numproj = K,
        numperm = 0
      )$stat
    }
    return(c(Tn_b >= Tn))
  }
  RES.cvm <- c(RES.cvm, mean((rowMeans(res) < alpha)))
  cat('\nDone~~~\n')
}

stopCluster(cl)
t2 <- Sys.time()
print((t2-t1))
###################### END #############################

plot(d.choose, RES.cvm, type = 'l')


save.image('Sim_d_res.RData')






