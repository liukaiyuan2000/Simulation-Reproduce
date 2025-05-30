rm(list = ls())
## some settings ----
library(splines)
library(MASS)
library(purrr)
library(doSNOW)
library(tcltk)
nsim = 10
burn_in = 500
mc_num = 1500
p = 2
K = 4
mu.c = 5
mu.r = 5
gamma.c = c(-2, 0.2)
gamma.r = c(-2, -0.2)
beta.true = c(0.5, 1)

progress <- function(n){
  pb <- txtProgressBar(max = nsim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)
######

## Functions ----
#' @description
#' Function for data generate process
#'
#' @param n sample size
#' @param cases simulation cases, 'additive' and 'linear'
#' @export x predictive variables
#' @export y response variable
#' @export z covariate about variance
#' @export B B-spline function matrix
DGP <- function(n, cases){
  if(cases == 'additive'){
    x1.c <- runif(n, -3, 3)
    x2.c <- runif(n)
    x.c  <- cbind(x1.c, x2.c)
    z.c  <- cbind(1, runif(n, -1, 1))

    sig2.c <- drop(exp(z.c %*% gamma.c))
    eps.c  <- rnorm(n) * sqrt(sig2.c)
    y.c    <- mu.c + x1.c + sin(2*pi*x2.c) + eps.c
    B.c    <- list()
    for(i in 1:p){
      B.c[[i]] <- scale(bs(x.c[, i], degree = 3, knots = mean(x.c[, i])), scale = F)
    }

    x1.r <- runif(n)
    x2.r <- runif(n)
    x.r  <- cbind(x1.r, x2.r)
    z.r  <- cbind(1, runif(n, 1, 2))

    sig2.r <- drop(exp(z.r %*% gamma.r))
    eps.r  <- rnorm(n) * sqrt(sig2.r)
    y.r    <- mu.r + cos(2*pi*x1.r) + sin(2*pi*x2.r) + eps.r
    B.r    <- list()
    for(i in 1:p){
      B.r[[i]] <- scale(bs(x.r[, i], degree = 3, knots = mean(x.r[, i])), scale = F)
    }
  } else if(cases == 'linear'){
    x1.c <- runif(n, -3, 3)
    x2.c <- runif(n, -3, 3)
    x.c  <- cbind(x1.c, x2.c)
    z.c  <- cbind(1, runif(n, -1, 1))

    sig2.c <- drop(exp(z.c %*% gamma.c))
    eps.c  <- rnorm(n) * sqrt(sig2.c)
    y.c    <- mu.c + drop(x.c %*% beta.true) + eps.c
    B.c    <- list()
    for(i in 1:p){
      B.c[[i]] <- scale(bs(x.c[, i], degree = 3, knots = mean(x.c[, i])), scale = F)
    }

    x1.r <- runif(n)
    x2.r <- runif(n)
    x.r  <- cbind(x1.r, x2.r)
    z.r  <- cbind(1, runif(n, 1, 2))

    sig2.r <- drop(exp(z.r %*% gamma.r))
    eps.r  <- rnorm(n) * sqrt(sig2.r)
    y.r    <- mu.r + cos(2*pi*x1.r) + sin(2*pi*x2.r) + eps.r
    B.r    <- list()
    for(i in 1:p){
      B.r[[i]] <- scale(bs(x.r[, i], degree = 3, knots = mean(x.r[, i])), scale = F)
    }
  }

  return(
    list(
      x.c = x.c, y.c = y.c,
      z.c = z.c, B.c = B.c,
      x.r = x.r, y.r = y.r,
      z.r = z.r, B.r = B.r
    )
  )
}

#' @description
#' Function for likelihood
#'
#' @param gamma
#' @param alpha
#' @param mu
#' @param x
#' @param y
#' @param z
#' @param B
#' @param gamma0 initial value of gamma
#' @export - likelihood value
p_gamma <- function(
    gamma, alpha, mu, x, y, z, B,
    gamma0 = c(0, 0), B_gamma_inv = diag(2)
){
  n = length(y)
  Sig <- diag(drop(exp(z %*% gamma)))
  temp1 <- rep(0, n)
  for(i in 1:p){
    temp1 <- temp1 + drop(B[[i]] %*% alpha[[i]])
  }
  f <- det(Sig)^(-1/2)*
    exp(
      -1/2*(
        t(y - temp1 - mu) %*% solve(Sig) %*% (y - temp1 - mu) +
          t(gamma - gamma0) %*% B_gamma_inv %*% (gamma - gamma0)
      )
    )
  return(drop(f))
}
#' @description
#' Function for Gibbs sampling procedure
#'
#' @param x
#' @param y
#' @param z
#' @param B
#' @param burn_in the number of burn_in period
#' @param mc_num the number of Markov chains
#' @param gamma0 initial value of gamma
#' @param mu0 initial value of mu
#' @param B_gamma
#' @param B_gamma_inv inverse of B_gamma
#' @export - a list contains the results for the parameters \mu, \gamma, and \alpha
Gibbs_sample <- function(
    x, z, y, B, burn_in, mc_num,
    gamma0 = c(0, 0), mu0 = 0,
    B_gamma = diag(2), B_gamma_inv = diag(2)
){
  alpha0 = list()
  for(i in 1:p){
    alpha0[[i]] = rep(0, K)
  }
  a_tau = b_tau = 1

  N = burn_in + mc_num

  gamma.old = gamma0
  mu.old    = mu0
  Sigma.old = diag(drop(exp(z %*% gamma.old)))

  alpha.iter    <- alpha0
  gamma.process <- matrix(0, N, 2)
  mu.process    <- rep(0, N)
  alpha.process <- list()
  for(i in 1:p){
    alpha.process[[i]] <- matrix(0, N, K)
  }
  p_acc <- rep(0, N)

  for(s in 1:N){
    cat(s, '\r')
    Sigma.old_inv <- solve(Sigma.old)
    # calculate tau2.new
    tau2.new <- list()
    for(i in 1:p){
      tau2.new[[i]] <- 1/rgamma(1, K/2 + a_tau, 1/2*(t(alpha.iter[[i]] - alpha0[[i]]) %*% (alpha.iter[[i]] - alpha0[[i]]) + 2*b_tau))
    }

    # calculate alpha.new
    y_minus_mu.old <- y - mu.old
    for(i in 1:p){
      temp1 <- rep(0, n)
      for(j in 1:p){
        temp1 <- temp1 + drop(B[[j]] %*% alpha.iter[[j]]) * (j != i)
      }
      b_alpha.star    <- solve(1/tau2.new[[i]]*diag(K) + t(B[[i]]) %*% Sigma.old_inv %*% B[[i]])
      alpha0.star     <- drop(b_alpha.star %*% (1/tau2.new[[i]]*diag(K) %*% alpha0[[i]] + t(B[[i]]) %*% Sigma.old_inv %*% (y_minus_mu.old - temp1)))
      alpha.iter[[i]] <- mvrnorm(1, alpha0.star, b_alpha.star)
      alpha.process[[i]][s, ] <- alpha.iter[[i]]
    }

    # calculate mu.new
    temp2 <- rep(0, n)
    for(j in 1:p){
      temp2 <- temp2 + drop(B[[j]] %*% alpha.iter[[j]])
    }
    mu.sig <- 1 / sum(diag(Sigma.old_inv))
    mu.mu  <- diag(Sigma.old_inv) %*% (y - temp2) * mu.sig
    mu.new <- rnorm(1, mu.mu, sqrt(mu.sig))

    # calculate gamma.new, MH Algorithm
    Omega_gamma <- 1/2 * t(z) %*% diag(drop((y - mu.new - temp2)^2 / exp(z %*% gamma.old))) %*% z + B_gamma_inv
    sig_gamma2  <- 1.5
    mc.new <- mvrnorm(1, gamma.old, sig_gamma2*solve(Omega_gamma))
    p_acc[s]  <- min(
      1, p_gamma(mc.new, alpha.iter, mu.new, x, y, z, B) /
        p_gamma(gamma.old, alpha.iter, mu.new, x, y, z, B)
    )
    gamma.new <- gamma.old + (mc.new - gamma.old)*(runif(1) < p_acc[s])
    gamma.process[s, ] <- gamma.old
    mu.process[s]      <- mu.old

    gamma.old = gamma.new
    Sigma.old = diag(drop(exp(z %*% gamma.new)))
    mu.old    = mu.new
  }
  return(
    list(
      gamma.process = gamma.process,
      alpha.process = alpha.process,
      mu.process    = mu.process
    )
  )
}

#' @description
#' Function of RMSEH
#'
#' @param L_hat estimate of lower bound
#' @param L true value of lower bound
#' @param U_hat estimate of upper bound
#' @param U true value of upper bound
#' @export - the value of RMSEH
RMSEH = function(L_hat, L, U_hat, U){
  temp = (sqrt(mean((L_hat- L)^2) + mean((U_hat- U)^2)))/2
  return(temp)
}

#' @description
#' Function of MAE
#'
#' @param L_hat estimate of lower bound
#' @param L true value of lower bound
#' @param U_hat estimate of upper bound
#' @param U true value of upper bound
#' @export - the value of MAE
MAE = function(L_hat, L, U_hat, U){
  temp = abs(L_hat - L) + abs(U_hat - U)
  return(mean(temp)/2)
}

#' @description
#' Function of AR
#'
#' @param L_hat estimate of lower bound
#' @param L true value of lower bound
#' @param U_hat estimate of upper bound
#' @param U true value of upper bound
#' @export - the value of AR
AR = function(L_hat, L, U_hat, U){
  temp <- (pmin(U_hat, U) - pmax(L_hat, L)) / (pmax(U_hat, U) - pmin(L_hat, L))
  return(mean(temp))
}

#' @description
#' Function for bind (only works for parallel computing)
#'
#' @param ... results (unspecified number)
#' @export - a list contains the results for the parameters \mu, \gamma, and \alpha
bind_fun<-function(...){
  alpha.store.c <- lapply(list(...), function(x) x$alpha.res.c)
  result.c <- NULL
  for (i in 1:length(alpha.store.c)) {
    if (is.null(result.c)) {
      result.c <- alpha.store.c[[i]]
    } else {
      result.c <- Map(rbind, result.c, alpha.store.c[[i]])
    }
  }
  alpha.store.r <- lapply(list(...), function(x) x$alpha.res.r)
  result.r <- NULL
  for (i in 1:length(alpha.store.r)) {
    if (is.null(result.r)) {
      result.r <- alpha.store.r[[i]]
    } else {
      result.r <- Map(rbind, result.r, alpha.store.r[[i]])
    }
  }
  data <- list(
    mu.res.c = do.call(c, lapply(list(...), function(x) x$mu.res.c)),
    gamma.res.c = do.call(rbind, lapply(list(...), function(x) x$gamma.res.c)),
    alpha.res.c = result.c,
    mu.res.r = do.call(c, lapply(list(...), function(x) x$mu.res.r)),
    gamma.res.r = do.call(rbind, lapply(list(...), function(x) x$gamma.res.r)),
    alpha.res.r = result.r,
    rmseh.res = do.call(c, lapply(list(...), function(x) x$rmseh.res)),
    mae.res = do.call(c, lapply(list(...), function(x) x$mae.res)),
    ar.res = do.call(c, lapply(list(...), function(x) x$ar.res))
  )
  return(data)
}
#####

