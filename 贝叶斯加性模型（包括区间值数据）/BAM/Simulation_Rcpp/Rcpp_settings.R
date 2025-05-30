rm(list = ls())
## some settings ----
library(splines)
library(MASS)
library(purrr)
library(Rcpp)
library(doSNOW)
library(tcltk)
library(rBAM)
nsim = 50
burn_in = 500
mc_num = 1500
p = 3
mu = 1
gamma = c(-0.5, 0.5)
gamma0 = c(-0.5, 0.5)
mu0 = 0
K = 4
alpha0 = list()
for(i in 1:p){
  alpha0[[i]] = rep(0, K)
}
B_gamma = diag(2)
B_gamma_inv = solve(B_gamma)
a_tau = b_tau = 1

progress <- function(n){
  pb <- txtProgressBar(max = nsim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)
######

#' @description
#' Function for bind (only works for parallel computing)
#'
#' @param ... results (unspecified number)
#' @export - a list contains the results for the parameters \mu, \gamma, and \alpha
bind_fun<-function(...){
  alpha.store <- lapply(list(...), function(x) x$alpha.res)
  result <- NULL
  for (i in 1:length(alpha.store)) {
    if (is.null(result)) {
      result <- alpha.store[[i]]
    } else {
      result <- Map(rbind, result, alpha.store[[i]])
    }
  }
  data <- list(
    mu.res = do.call(c, lapply(list(...), function(x) x$mu.res)),
    gamma.res = do.call(rbind, lapply(list(...), function(x) x$gamma.res)),
    alpha.res = result
  )
  return(data)
}


