rm(list = ls())
library(doSNOW)
library(tcltk)
library(Sinica2022)
no_cores <- 12
progress <- function(n){
  pb <- txtProgressBar(max = nsim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)

m = 25
alpha = 0.05
nsim = 500
rho = 300

Sim.fun <- function(n, C, alpha, m, setting, nsim, rho){
  result <- foreach(
    sim = 1:nsim, .options.snow = opts, .packages = 'Sinica2022',
    .combine = 'rbind', .errorhandling = "remove"
  )  %dopar% {
    ## data and kernel parameters
    data = DGP(n, C, setting)
    X.tilde = data$X.tilde
    Y.tilde = data$Y.tilde
    Z = data$Z
    U = data$U
    u.mat = matrix(rep(U, n), n, n, byrow = F) - matrix(rep(U, n), n, n, byrow = T)
    h = 2.34 * sd(U) * n^(-1/3)
    kern.u = 0.75 * (1 - (u.mat/h)^2) * (abs(u.mat) <= h) / h
    ## calibrating
    X.hat = apply(X.tilde, 2, \(x) cali.func(x, u.mat, kern.u))
    Y.hat = cali.func(Y.tilde, u.mat, kern.u)
    ## LS estimate and residual
    beta_e.hat = beta.hat.func(X.hat, Z, Y.hat, setting)
    beta.hat = beta_e.hat$beta.hat
    e.hat = beta_e.hat$e.hat
    V.hat = cbind(X.hat, Z)
    J = ncol(V.hat)
    theta.uniform <- matrix(runif(m*J, -1, 1), J, m)
    theta.normal <- matrix(rnorm(m*J), J, m)
    ## test statistics
    Tcvm.res  <- Tcvm(X.hat, Z, e.hat)
    Tcvm.hat  <- Tcvm.res$Tn
    Tcvmu.hat <- Tcvm.rand(V.hat, theta.uniform, e.hat)
    Tcvmn.hat <- Tcvm.rand(V.hat, theta.normal, e.hat)
    Tksu.hat  <- Tks_rand_rcpp(V.hat, theta.uniform, e.hat)
    Tksn.hat  <- Tks_rand_rcpp(V.hat, theta.normal, e.hat)
    ## bootstrap
    A.hat <- Tcvm.res$A
    e.star = bootstrap.res(X.hat, Z, Y.hat, e.hat, setting, rho)
    Tcvm.star  <- diag(t(e.star) %*% A.hat %*% e.star) / n^2
    Tcvmu.star <- Tcvm.rand(V.hat, theta.uniform, e.star)
    Tcvmn.star <- Tcvm.rand(V.hat, theta.normal, e.star)
    Tksu.star  <- Tks_rand_rcpp(V.hat, theta.uniform, e.star)
    Tksn.star  <- Tks_rand_rcpp(V.hat, theta.normal, e.star)

    res1 = (Tcvm.hat >= sort(Tcvm.star)[(1 - alpha)*rho])
    res2 = (Tcvmu.hat >= sort(Tcvmu.star)[(1 - alpha)*rho])
    res3 = (Tcvmn.hat >= sort(Tcvmn.star)[(1 - alpha)*rho])
    res4 = (Tksu.hat >= sort(Tksu.star)[(1 - alpha)*rho])
    res5 = (Tksn.hat >= sort(Tksn.star)[(1 - alpha)*rho])

    return(c(res1, res2, res3, res4, res5))
  }
  RES <- colMeans(result, na.rm = T)
  return(RES)
}





