rm(list = ls())
## some settings ----
library(splines)
library(MASS)
library(purrr)
library(doSNOW)
library(tcltk)
nsim = 1
burn_in = 1000
mc_num = 2000
p = 3
mu = 10
gamma = c(-0.5, 1)
gamma0 = c(0, 0)
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
#' Function for data generate process
#' 
#' @param n sample size
#' @export x predictive variables
#' @export y response variable
#' @export z covariate about variance
#' @export B B-spline function matrix
DGP <- function(n){
  x1 <- runif(n, -1, 1)
  x2 <- runif(n)
  x3 <- runif(n)
  x  <- cbind(x1, x2, x3)
  z  <- cbind(1, runif(n, -1, 1))
  
  sig2 <- drop(exp(z %*% gamma))
  eps  <- rnorm(n) * sqrt(sig2)
  y    <- x1 + mu + sin(2*pi*x2) + cos(2*pi*x3) + eps
  B    <- list()
  for(i in 1:p){
    B[[i]] <- scale(bs(x[, i], degree = 3, knots = mean(x[, i])), scale = F)
  }
  return(
    list(
      x1 = x1, x2 = x2, x = x, 
      y = y, z = z, B = B, 
      sig2 = sig2, eps = eps
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
#' @export - likelihood value
p_gamma <- function(gamma, alpha, mu, x, y, z, B){
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
#' @param n sample size
#' @param burn_in the number of burn_in period
#' @param mc_num the number of Markov chains
#' @export - a list contains the results for the parameters \mu, \gamma, and \alpha
Gibbs_sample <- function(n, burn_in, mc_num){
  data <- DGP(n)
  x <- data$x
  z <- data$z
  y <- data$y
  B <- data$B
  
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
      data = data, 
      gamma.process = gamma.process, 
      alpha.process = alpha.process, 
      mu.process    = mu.process,
      y = y, 
      B = B
    )
  )
}
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