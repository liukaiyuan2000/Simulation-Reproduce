# A Kernel Method for the Two-Sample-Problem
rm(list = ls())
set.seed(123)
library(PMMD)
library(doSNOW)
library(tcltk)
library(copula)
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

no_cores <- 10
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
        ## Cx=Cy, Fx!=Fy
        cop <- normalCopula(dim = d, param = 0.5)
        element.X <- list(mean = 0, sd = 1)
        element.Y <- list(df = 3)
        para.X <- vector("list", d)
        para.Y <- vector("list", d)
        for (i in 1:d) {
          para.X[[i]] <- element.X
          para.Y[[i]] <- element.Y
        }
        mv.X <- mvdc(cop, rep("norm", d), para.X)
        mv.Y <- mvdc(cop, rep("t", d), para.Y)
        X <- rMvdc(n, mv.X)
        Y <- rMvdc(n, mv.Y)
      },
      {
        ## Cx!=Cy, Fx=Fy
        cop.X <- normalCopula(dim = d, param = 0.8)
        cop.Y <- tCopula(dim = d, param = 0.2, df = 3)
        element <- list(mean = 0, sd = 1)
        para <- vector("list", d)
        for (i in 1:d) {
          para[[i]] <- element
        }
        mv.X <- mvdc(cop.X, rep("norm", d), para)
        mv.Y <- mvdc(cop.Y, rep("norm", d), para)
        X <- rMvdc(n, mv.X)
        Y <- rMvdc(n, mv.Y)
      },
      {
        ## Cx=Cy, Fx1=Fy1, Fx2!=Fy2
        cop <- normalCopula(dim = d, param = 0.5)
        element.X <- list(mean = 0, sd = 1)
        element.Y <- list(df = 3)
        para.X <- vector("list", d)
        para.Y <- vector("list", d)
        for (i in 1:d) {
          para.X[[i]] <- element.X
          if(i <= floor(d/2)){
            para.Y[[i]] <- element.Y
          } else{
            para.Y[[i]] <- element.X
          }
        }
        mv.X <- mvdc(cop, rep("norm", d), para.X)
        mv.Y <- mvdc(cop, c(rep("t", floor(d/2)), rep("norm", d - floor(d/2))), para.Y)
        X <- rMvdc(n, mv.X)
        Y <- rMvdc(n, mv.Y)
      }, 
      {
        ## Cx!=Cy, Fx1=Fy1, Fx2!=Fy2
        cop.X <- normalCopula(dim = d, param = 0.8)
        cop.Y <- tCopula(dim = d, param = 0.2, df = 3)
        element.X <- list(mean = 0, sd = 1)
        element.Y <- list(df = 3)
        para.X <- vector("list", d)
        para.Y <- vector("list", d)
        for (i in 1:d) {
          para.X[[i]] <- element.X
          if(i <= floor(d/2)){
            para.Y[[i]] <- element.Y
          } else{
            para.Y[[i]] <- element.X
          }
        }
        mv.X <- mvdc(cop.X, rep("norm", d), para.X)
        mv.Y <- mvdc(cop.Y, c(rep("t", floor(d/2)), rep("norm", d - floor(d/2))), para.Y)
        X <- rMvdc(n, mv.X)
        Y <- rMvdc(n, mv.Y)
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
    sim = 1:Sim, .options.snow = opts, 
    .packages = c('PMMD', 'copula'), .combine = 'rbind'
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

for(setting in 1:4){
  RES.cvm = c()
  RES.ks  = c()
  for(n in n.choose){
    for(d in d.choose){
      cat("Case ", setting, "; n = ", n, "; d = ", d, "\n", sep = "")
      res = Sim_fun(n, d, setting, B)
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
  RES.sum[[paste("setting = ", setting)]] <- list(
    CvM = res.cvm.dat, KS = res.ks.dat
  )
}

stopCluster(cl)
t2 <- Sys.time()
print((t2-t1))
###################### END #############################

RES.sum
save.image('copula_normal_t3_res2.RData')


RES.dat = rbind(
  cbind(RES.sum$`setting =  1`$CvM, RES.sum$`setting =  1`$KS), 
  cbind(RES.sum$`setting =  2`$CvM, RES.sum$`setting =  2`$KS), 
  cbind(RES.sum$`setting =  3`$CvM, RES.sum$`setting =  3`$KS), 
  cbind(RES.sum$`setting =  4`$CvM, RES.sum$`setting =  4`$KS)
)
rownames(RES.dat) <- NULL
RES.dat = as.matrix(cbind(" " = rep(c('$n=50$', '$n=100$', '$n=200$'), 4), RES.dat))
stargazer(cbind(NA, RES.dat), summary = F, digits = 3, table.placement = "H")

