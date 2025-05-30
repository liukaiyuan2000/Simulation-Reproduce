# A Kernel Method for the Two-Sample-Problem
rm(list = ls())
set.seed(123)
library(PMMD)
library(doSNOW)
library(tcltk)
library(MASS)
library(mvtnorm)
library(stargazer)
setwd('D:/Desktop/MMD')
############# SETTING AND FUNCTION ##################
n.choose <- c(50, 100, 200)
d.choose <- c(5, 10, 20, 40)

B = 200
Sim = 1000
alpha = 0.05
res_colname <- c("d=5", 'd=10', 'd=20', 'd=40')
res_rowname <- c('n=50', 'n=100', 'n=200')

no_cores <- 12
progress <- function(n){
  pb <- txtProgressBar(max = Sim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)
Sim_fun <- function(n, d, B = 200, K = 10, alpha = 0.05){
  m = n
  # Temp = abs(outer(1:d, 1:d, '-'))
  # Sig = (Temp <= 2) * 0.2 + 0.8*diag(d)
  Sig = 0.5^abs(outer(1:d, 1:d, '-'))
  result <- foreach(
    sim = 1:Sim, .options.snow = opts, 
    .packages = c('PMMD', 'MASS', 'mvtnorm'), .combine = 'rbind'
  )  %dopar% {
    
    X = rmvnorm(n, rep(0, d), Sig)
    Y = rmvt(n, Sig, df = 3)
    
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

RES.cvm = c()
RES.ks  = c()
for(n in n.choose){
  for(d in d.choose){
    cat("n = ", n, "; d = ", d, "\n", sep = "")
    res = Sim_fun(n, d)
    RES.cvm <- c(RES.cvm, res$CvM)
    RES.ks  <- c(RES.ks, res$KS)
    cat('\nDone~~~\n')
  }
}
res.cvm.dat <- as.data.frame(
  matrix(RES.cvm, length(n.choose), length(d.choose), byrow = T),
  row.names = res_rowname
)
res.ks.dat <- as.data.frame(
  matrix(RES.ks, length(n.choose), length(d.choose), byrow = T),
  row.names = res_rowname
)
colnames(res.cvm.dat) = res_colname
colnames(res.ks.dat) = res_colname
RES.sum <- list(
  CvM = res.cvm.dat, KS = res.ks.dat
)

stopCluster(cl)
t2 <- Sys.time()
print((t2-t1))
###################### END #############################

RES.sum
save.image('multinormalres3.RData')


RES.dat = cbind(RES.sum$CvM, RES.sum$KS)
rownames(RES.dat) <- NULL
RES.dat = as.matrix(cbind(" " = c('$n=50$', '$n=100$', '$n=200$'), RES.dat))
stargazer(RES.dat, summary = F, digits = 3, table.placement = "H")





