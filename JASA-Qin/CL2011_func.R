rm(list = ls())
library(MASS)
library(Matrix)

gen_sigma <- function(p, model){
  switch(
    model, 
    {
      p_half <- p / 2
      temp <- 1 - abs(outer(1:p_half, 1:p_half, '-')) / 10
      A1 <- temp * (temp >= 0)
      A2 <- diag(4, p_half)
      sigma <- as.matrix(bdiag(A1, A2))
    }, 
    {
      p_half <- p / 2
      B = matrix(rbinom(p_half^2, 1, 0.2)*runif(p_half^2, 0.3, 0.8), p_half, p_half)
      rowscols = which(upper.tri(B), arr.ind = TRUE)
      B = sparseMatrix(
        i = rowscols[,1], j = rowscols[,2], x = B[upper.tri(B)], 
        symmetric = TRUE, dims = c(p_half, p_half)
      )
      diag(B) = rbinom(p_half, 1, 0.2)*runif(p_half, 0.3, 0.8)
      lambda = eigen(B, only.values = T)$values
      if(is.complex(lambda)){
        eps = max(-min(Re(lambda)), 0) + 0.01
      } else {
        eps = max(-min(lambda), 0) + 0.01
      }
      A1 <- B + eps * diag(p_half)
      A2 <- diag(4, p_half)
      sigma <- as.matrix(bdiag(A1, A2))
    }
  )
  return(sigma)
}

cal_risk <- function(cov_1, cov_2, norm_type = "operator"){
  stopifnot(all(dim(cov_1) == dim(cov_2)))
  diff <- cov_1 - cov_2
  if (norm_type == "operator") {
    risk <- norm(diff, "2")
  } else if (norm_type == "l1") {
    risk <- norm(diff, "1")
  } else if (norm_type == "frobenius") {
    risk <- norm(diff, "F")
  } else {
    stop("norm_type must be 'spectral', 'l1' or 'frobenius'")
  }
  return(risk)
}

thresholding <- function(cov, threshold, s_type = "hard"){
  stopifnot(all(dim(cov) == dim(threshold)))
  if (s_type == "hard") {
    cov[abs(cov) < threshold] <- 0
  } else if (s_type == "adaptive lasso") {
    temp = 1 - abs(threshold / cov)^4
    idx = (temp > 0)
    cov = cov * temp * idx
  } else {
    stop("method must be 'adaptive lasso' or 'hard'")
  }
  
  return(cov)
}

cal_theta_hat <- function(X){
  n <- nrow(X)
  p <- ncol(X)
  sigma_hat <- cov(X)
  X_colmean <- colMeans(X)
  theta_hat <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    x_i_bar <- X_colmean[i]
    for (j in 1:p) {
      x_j_bar <- X_colmean[j]
      theta_hat[i, j] <- mean(
        ((X[, i] - x_i_bar)*(X[, j] - x_j_bar) - sigma_hat[i, j])^2
      )
    }
  }
  return(theta_hat)
}

thre_cv_func <- function(
    X, H = 5, method = "adaptive", s_type = "hard", delta.len = 10
){
  p <- ncol(X)
  n <- nrow(X)
  max_cov <- max(abs(cov(X)))
  delta_candidates <- seq(0, 4, length.out = delta.len)
  
  if (method == "universal") {
    threshold_candidates <- lapply(delta_candidates, function(delta) {
      matrix(delta * sqrt(log(p) / n), nrow = p, ncol = p)
    })
  } else if (method == "adaptive") {
    theta_hat <- cal_theta_hat(X)
    threshold_candidates <- lapply(delta_candidates, function(delta) {
      delta * sqrt(theta_hat * log(p) / n)
    })
  } else if (method == "adaptive2") {
    theta_hat <- cal_theta_hat(X)
    return(2 * sqrt(theta_hat * log(p) / n))
  } else {
    stop("method must be 'universal', 'adaptive' or 'adaptive2'")
  }
  
  risks <- numeric(length(threshold_candidates))
  ## index for CV
  folds <- sample(rep(1:H, each = n / H))
  
  for (idx in seq_along(threshold_candidates)) {
    threshold <- threshold_candidates[[idx]]
    risk <- 0
    for (fold in 1:H) {
      train_idx <- which(folds != fold)
      test_idx <- which(folds == fold)
      cov_n_1 <- cov(X[train_idx, ])
      cov_n_2 <- cov(X[test_idx, ])
      cov_n_1 <- thresholding(cov_n_1, threshold, s_type)
      risk <- risk + cal_risk(cov_n_1, cov_n_2, norm_type = "frobenius")^2
    }
    risks[idx] <- risk / H
  }
  
  best_threshold <- threshold_candidates[[which.min(risks)]]
  return(best_threshold)
}

cal_F_statistics <- function(X, Y) {
  n <- nrow(X)
  k <- length(unique(Y))
  x_bar <- colMeans(X)
  
  numerator <- 0
  for (m in unique(Y)) {
    df_m <- X[Y == m, ]
    x_bar_m <- colMeans(df_m)
    n_m <- nrow(df_m)
    numerator <- numerator + n_m * (x_bar_m - x_bar)^2
  }
  numerator <- numerator / (k - 1)
  
  denominator <- 0
  for (m in unique(Y)) {
    df_m <- X[Y == m, ]
    sigma_m_square <- apply(df_m, 2, var)
    n_m <- nrow(df_m)
    denominator <- denominator + (n_m - 1) * sigma_m_square
  }
  denominator <- denominator / (n - k)
  
  res <- numerator / denominator
  return(res)
}
