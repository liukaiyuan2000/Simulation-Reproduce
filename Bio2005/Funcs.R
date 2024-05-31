rm(list = ls())
library(purrr)
library(doSNOW)
library(tcltk)
library(Bio2005)
nsim = 1000
n = m = 2^(2:8)
idx <- length(m)
p_t = p_w = matrix(0, length(n), length(m))
no_cores <- 5
progress <- function(n){
  pb <- txtProgressBar(max = nsim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)

Sim_fun <- function(case, m, n){
  m_long <- rep(m, each = idx)
  n_long <- rep(n, idx)
  rho <- switch(case, 0, 0.1)
  result <- foreach(
    sim = 1:nsim, .options.snow = opts, .packages = c('purrr', 'Bio2005')
  )  %dopar% {
    N_long = n_long + 1
    y <- map2(N_long, m_long, \(x, y) DGP(x, y, rho))
    
    p_t1 <- map_dbl(y, \(x) Tnm_rcpp(x)[2])
    p_w1 <- map_dbl(y, \(x) Wnm_rcpp(x)[2])
    return(
      data.frame(p_t1 = p_t1, p_w1 = p_w1)
    )
  }
  return(result)
}

Res_fun <- function(RES_sum){
  T1 <- do.call(rbind, map(1:nsim, \(x) RES_sum[[x]]$p_t1))
  T2 <- do.call(rbind, map(1:nsim, \(x) RES_sum[[x]]$p_w1))
  p_t <- matrix(apply(T1 < 0.05, 2, mean), idx, idx, byrow = T)
  p_w <- matrix(apply(T2 < 0.05, 2, mean), idx, idx, byrow = T)
  p_w[!upper.tri(p_w, diag = TRUE)] <- NA
  colnames(p_t) <- paste("n =",n)
  rownames(p_t) <- paste("m =",m)
  colnames(p_w) <- paste("n =",n)
  rownames(p_w) <- paste("m =",m)
  return(
    list(p_t = p_t, p_w = p_w)
  )
}


