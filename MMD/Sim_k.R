# A Kernel Method for the Two-Sample-Problem
rm(list = ls())
set.seed(123)
library(PMMD)
library(doSNOW)
library(tcltk)
library(stargazer)
library(MASS)
setwd('D:/Desktop/MMD')
############# SETTING AND FUNCTION ##################
n = 40
d.choose <- c(5, 10, 20, 40)
K.choose <- c(10, 50, 100, 150, 200)

B = 200
Sim = 1000
alpha = 0.05
res_rowname <- c("$d=5$", '$d=10$', '$d=20$', '$d=40$')
res_colname <- c('$K=10$', '$K=50$', '$K=100$', '$K=150$','$K=200$')

no_cores <- 12
progress <- function(n){
  pb <- txtProgressBar(max = Sim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)
Sim_fun <- function(n, d, setting, B = 100, K = 10, alpha = 0.05){
  m = n
  DGP <- function(n, d, setting){
    m = n
    switch(
      setting,
      {
        X <- matrix(rnorm(n*d), n, d)
        Y <- matrix(rnorm(m*d), m, d)
      },
      {
        X <- matrix(rnorm(n*d), n, d)
        Y <- matrix(rnorm(m*d, 0.1), n, d)
      },
      {
        X <- matrix(rnorm(n*d), n, d)
        Y <- mvrnorm(m, rep(0, d), matrix(0.1, d, d) + 0.9*diag(d))
      }
    )
    return(
      list(
        X = X,
        Y = Y
      )
    )
  }
  result <- foreach(
    sim = 1:Sim, .options.snow = opts, .packages = c('PMMD', 'MASS'), .combine = 'rbind'
  )  %dopar% {
    dat = DGP(n, d, setting)
    X = dat$X
    Y = dat$Y
    
    Tn <- meammd(
      X = X, Y = Y, pval = TRUE, type = "proj", numproj = K,
      numperm = 0
    )$stat
    KSTn <- sqrt(maxmmd(
      X = X, Y = Y, pval = TRUE, type = "proj", numproj = K,
      numperm = 0
    ))
    Tn_b = rep(0, B)
    KSTn_b = rep(0, B)
    for (b in 1:B) {
      Z <- rbind(X, Y)
      set <- sample(1:(n+m), n)
      
      X1 <- Z[set, ]
      Y1 <- Z[-set,]
      
      Tn_b[b] <- meammd(
        X = X1, Y = Y1, pval = TRUE, type = "proj", numproj = K,
        numperm = 0
      )$stat
      KSTn_b[b] <- sqrt(maxmmd(
        X = X1, Y = Y1, pval = TRUE, type = "proj", numproj = K,
        numperm = 0
      ))
    }
    return(c(Tn_b >= Tn, KSTn_b >= KSTn))
  }
  return(
    data.frame(
      CvM = mean((rowMeans(result[, 1:B]) < alpha)),
      KS  = mean((rowMeans(result[, (B+1):(2*B)]) < alpha))
    )
  )
}
#######################################################

##################### Simulation ######################
t1 <- Sys.time()
cl <- makeSOCKcluster(no_cores)
registerDoSNOW(cl)
RES.sum = list()

for(setting in 1:3){
  RES.cvm = c()
  RES.ks  = c()
  for(d in d.choose){
    for(K in K.choose){
      cat("Case ", setting, "; K = ", K, "; d = ", d, "\n", sep = "")
      res = Sim_fun(n, d, setting, B, K)
      RES.cvm <- c(RES.cvm, res$CvM)
      RES.ks  <- c(RES.ks, res$KS)
      cat('\nDone~~~\n')
    }
  }
  
  res.cvm.dat <- as.data.frame(
    matrix(RES.cvm, length(d.choose), length(K.choose), byrow = T),
    row.names = res_rowname
  )
  res.ks.dat <- as.data.frame(
    matrix(RES.ks, length(d.choose), length(K.choose), byrow = T),
    row.names = res_rowname
  )
  colnames(res.cvm.dat) = res_colname
  colnames(res.ks.dat) = res_colname
  RES.sum[[paste("setting = ", setting)]] <- list(
    CvM = res.cvm.dat, KS = res.ks.dat
  )
}

stopCluster(cl)
t2 <- Sys.time()
print((t2-t1))
###################### END #############################

RES.sum

write.csv(RES.sum, 'K_res.csv')


RES.dat = rbind(
  cbind(RES.sum$`setting =  1`$CvM, RES.sum$`setting =  1`$KS), 
  cbind(RES.sum$`setting =  2`$CvM, RES.sum$`setting =  2`$KS), 
  cbind(RES.sum$`setting =  3`$CvM, RES.sum$`setting =  3`$KS)
)
rownames(RES.dat) <- NULL
RES.dat = as.matrix(cbind(" " = rep(c("$d=5$", '$d=10$', '$d=20$', '$d=40$'), 3), RES.dat))
stargazer(cbind(NA, RES.dat), summary = F, digits = 3, table.placement = "H")

