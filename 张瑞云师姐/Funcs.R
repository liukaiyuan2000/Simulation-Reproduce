rm(list = ls())
library(purrr)
library(doSNOW)
library(tcltk)
library(PFLSAM2024)
nsim = 1000
m <- 5
R <- 30
n = m * R
rho = 0.25
alpha = 0.05
B = 500

no_cores <- 6
progress <- function(n){
  pb <- txtProgressBar(max = nsim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)

Sim_cvm <- function(R, m, rho, case.DGP, case.var, alpha, B){
  result <- foreach(
    sim = 1:nsim, .options.snow = opts, .packages = c('purrr', 'PFLSAM2024'), .combine = 'c'
  ) %dopar% {
    n = m * R
    data <- DGP(R, m, rho, case.DGP, case.var)
    Wn <- data$Wn
    Zn <- data$Zn
    Xnt <- data$Xnt
    Yn <- data$Yn
    vec <- phi.hat(Xnt)
    Phi <- vec$Phi
    h_2_tilde <- vec$h_2_tilde
    paras <- para_hat_rcpp(Zn, Wn, Yn, Phi)
    theta.hat <- paras$theta_hat
    alpha.hat <- paras$alpha_hat
    test.val <- CTn_cvm_rcpp(Zn, Wn, Yn, Phi, Xnt, theta.hat, alpha.hat, h_2_tilde)
    varepsilon.hat <- test.val$varepsilon_hat
    CTn <- test.val$CT_n
    CTn.star <- bootstrap_cvm_rcpp(Zn, Wn, Phi, Xnt, theta.hat, alpha.hat, h_2_tilde, varepsilon.hat, B)
    return(1 * (mean(CTn.star > CTn) < alpha))
  }

  return(result)
}

Sim_ks <- function(R, m, rho, case.DGP, case.var, alpha, B){
  result <- foreach(
    sim = 1:nsim, .options.snow = opts, .packages = c('purrr', 'PFLSAM2024'), .combine = 'c'
  ) %dopar% {
    n = m * R
    data <- DGP(R, m, rho, case.DGP, case.var)
    Wn <- data$Wn
    Zn <- data$Zn
    Xnt <- data$Xnt
    Yn <- data$Yn
    vec <- phi.hat(Xnt)
    Phi <- vec$Phi
    h_2_tilde <- vec$h_2_tilde
    paras <- para_hat_rcpp(Zn, Wn, Yn, Phi)
    theta.hat <- paras$theta_hat
    alpha.hat <- paras$alpha_hat
    test.val <- CTn_ks_rcpp(Zn, Wn, Yn, Phi, Xnt, theta.hat, alpha.hat, h_2_tilde)
    varepsilon.hat <- test.val$varepsilon_hat
    CTn <- test.val$CT_n
    CTn.star <- bootstrap_ks_rcpp(Zn, Wn, Phi, Xnt, theta.hat, alpha.hat, h_2_tilde, varepsilon.hat, B)
    return(1 * (mean(CTn.star > CTn) < alpha))
  }

  return(result)
}





