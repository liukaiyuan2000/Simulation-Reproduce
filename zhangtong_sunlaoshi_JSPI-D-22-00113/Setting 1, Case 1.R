rm(list = ls())
Rcpp::sourceCpp("D:/Desktop/funcs.cpp")
#generate data function
DGP = function(n, cov.u, delta) {  
  const = (c(1:50) - 0.5) * pi  
  g = do.call(cbind, lapply(const^(-2), \(x) rnorm(n, 0, sqrt(x))))  
  phi = sqrt(2) * sin(const %*% t(tm))
  x = g %*% phi
  z1 = rnorm(n, -1, 1)
  z2 = runif(n, 0, pi)
  z = cbind(z1, z2)
  u = matrix(rnorm(n*2, 0, sqrt(cov.u)), n, 2)
  w = z + u
  y = rep(0, n)
  error = rnorm(n, 0, 0.5)
  for (i in c(1:n)) {  
    y[i] = g[i, 1] + 3 * g[i, 2] + z[i, ] %*% beta + error[i] + delta * sqrt(inpro(x[i, ], x[i, ], tm)) * abs(z[i, 1] - z[i, 2]) * 1.5  
  }
  return(list(Y = y, X = x, Z = z, W = w))
}


testlevel = 0.95
tm = c(1:100)*0.01 - 0.005
beta = c(-1.5, 1)
B = 300
m = 2
nsim = 500
RES = c(0, seq(1, 2.5, 0.5))

t1 <- Sys.time()
for(n in c(100, 200)){
  BE.res.store = c()
  NA.res.store = c()
  for(delta in c(0, seq(1, 2.5, 0.5))){
    BE.result = rep(0, nsim)
    NA.result = rep(0, nsim)
    for(sim in c(1:nsim)){
      cat(n, delta, sim, '\r')
      data = DGP(n, 0.8, delta)
      Y = data$Y
      Z = data$Z
      X = data$X
      W = data$W
      # Benchmark method
      BE.res = BE_stat(X, Y, Z, 2, tm)
      BE.stat = BE.res$stat
      BE.alphahatf = BE.res$alphahatf
      BE.betahat = BE.res$betahat
      BE.H = BE.res$H
      BE.S = BE.res$S
      BE.eigf = BE.res$eigf
      BE.ehatcun = BE.res$ehat
      BE.Kz = BE.res$Kz
      BE.bvalue = BE_boot(X, Z, BE.alphahatf, BE.betahat, BE.H, BE.S, BE.eigf, BE.ehatcun, BE.Kz, B, tm)
      # Naive method
      NA.res = BE_stat(X, Y, W, 2, tm)
      NA.stat = NA.res$stat
      NA.alphahatf = NA.res$alphahatf
      NA.betahat = NA.res$betahat
      NA.H = NA.res$H
      NA.S = NA.res$S
      NA.eigf = NA.res$eigf
      NA.ehatcun = NA.res$ehat
      NA.Kz = NA.res$Kz
      NA.bvalue = BE_boot(X, W, NA.alphahatf, NA.betahat, NA.H, NA.S, NA.eigf, NA.ehatcun, NA.Kz, B, tm)
      
      BE.result[sim] = (abs(BE.stat) > quantile(abs(BE.bvalue), testlevel))
      NA.result[sim] = (abs(NA.stat) > quantile(abs(NA.bvalue), testlevel))
    }
    BE.res.store = c(BE.res.store, mean(BE.result))
    NA.res.store = c(NA.res.store, mean(NA.result))
    # print(c("Benckmark", n, delta, mean(BE.result), testlevel))
    # print(c("Naive", n, delta, mean(NA.result), testlevel))
  }
  RES = rbind(RES, BE.res.store, NA.res.store)
}
t2 <- Sys.time()
print(t2-t1)

RES






