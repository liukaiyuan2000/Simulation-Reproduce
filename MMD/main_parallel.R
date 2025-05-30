#A Kernel Method for the Two-Sample-Problem##d.K都是单个写的
rm(list = ls())
# set.seed(123)
library(PMMD)
library(doSNOW)
library(tcltk)
############# SETTING AND FUNCTION ##################
n.choose <- c(40)
d.choose <- c(5, 10)

B = 100
Sim = 1000
alpha = 0.05
res_colname <- c("d = 5", 'd = 10')
res_rowname <- c('n = 40')

no_cores <- 15
progress <- function(n){
  pb <- txtProgressBar(max = Sim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)
DGP <- function(n, d, setting){
  m = n
  switch(
    setting,
    {
      X <- matrix(rnorm(n*d), n, d)
      Y <- rbind(matrix(rnorm(n*d/2, -1, 0.2), n/2, d), matrix(rnorm(n*d/2, 1, 0.2), n/2, d))
    },
    {
      X <- matrix(rnorm(n*d), n, d)
      Y <- matrix(rt(m*d, df = 3), m, d)
    },
    {
      X <- matrix(rnorm(n*d), n, d)
      Y <- matrix(rexp(m*d) - 1, m, d)
    }
  )
  return(
    list(
      X = X,
      Y = Y
    )
  )
}
Sim_fun <- function(n, d, setting, B = 100, K = 10, alpha = 0.05){
  m = n
  DGP <- function(n, d, setting){
    m = n
    switch(
      setting,
      {
        X <- matrix(runif(n*d, 0, 1), n, d)
        Y <- matrix(runif(m*d, 0, 0.98), m, d)
      },
      {
        X <- matrix(rnorm(n*d), n, d)
        Y <- matrix(rexp(m*d) - 1, m, d)
      },
      {
        X <- matrix(rnorm(n*d, 1, sqrt(2)), n, d)
        Y <- matrix(rchisq(m*d, 1), m, d)
      },
      {
        X <- matrix(rnorm(n*d), n, d)
        Y <- rbind(matrix(rnorm(n*d/2, -1, 0.2), n/2, d), matrix(rnorm(n*d/2, 1, 0.2), n/2, d))
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
    sim = 1:Sim, .options.snow = opts, .packages = 'PMMD', .combine = 'rbind'
  )  %dopar% {
    # X <- matrix(rnorm(n*d, 0, 1), n, d)
    # # Y <- matrix(rnorm(m*d, delta, 1 - rho), m, d)
    # Y <- matrix(rt(m*d, df = 5), m, d)

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
      CvM = mean((rowMeans(result) < alpha)[1:B]),
      KS  = mean((rowMeans(result) < alpha)[(B+1):(2*B)])
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
  for(n in n.choose){
    for(d in d.choose){
      cat("n = ", n, "; d = ", d, "\n", sep = "")
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




