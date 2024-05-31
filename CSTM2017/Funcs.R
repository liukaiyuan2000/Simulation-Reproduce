rm(list = ls())
library(purrr)
library(doSNOW)
library(tcltk)
library(CSTM2017)
nsim = 1000
N = 2^(2:7)
m = 2^(2:7)
idxm <- length(m)
idxN <- length(N)
p_t = p_w = matrix(0, length(N), length(m))
no_cores <- 12
progress <- function(n){
  pb <- txtProgressBar(max = nsim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)

Sim_fun <- function(case, m, N){
  m_long <- rep(m, each = idxN)
  N_long <- rep(N, idxm)
  result <- foreach(
    sim = 1:nsim, .options.snow = opts, .packages = c('purrr', 'CSTM2017')
  )  %dopar% {
    y <- map2(N_long, m_long, \(x, y) DGP(x, y, case))
    
    p_t1 <- map_dbl(y, \(x) Tnm_rcpp(x)[2])
    p_s1 <- map_dbl(y, \(x) Snm_rcpp(x)[2])
    p_trho1 <- map_dbl(y, \(x) Trho_rcpp(x)[2])
    p_t31 <- map_dbl(y, \(x) T3_rcpp(x)[2])
    return(
      data.frame(p_t1 = p_t1, p_s1 = p_s1, p_trho1 = p_trho1, p_t31 = p_t31)
    )
  }
  return(result)
}

Res_fun <- function(RES_sum){
  T1 <- do.call(rbind, map(1:nsim, \(x) RES_sum[[x]]$p_t1))
  T2 <- do.call(rbind, map(1:nsim, \(x) RES_sum[[x]]$p_s1))
  T3 <- do.call(rbind, map(1:nsim, \(x) RES_sum[[x]]$p_trho1))
  T4 <- do.call(rbind, map(1:nsim, \(x) RES_sum[[x]]$p_t31))
  p_t <- matrix(apply(T1 < 0.05, 2, mean), idxN, idxm, byrow = F)
  p_s <- matrix(apply(T2 < 0.05, 2, mean), idxN, idxm, byrow = F)
  p_trho <- matrix(apply(T3 < 0.05, 2, mean), idxN, idxm, byrow = F)
  p_t3 <- matrix(apply(T4 < 0.05, 2, mean), idxN, idxm, byrow = F)
  rownames(p_t) <- paste("N =", N)
  colnames(p_t) <- paste("m =", m)
  rownames(p_s) <- paste("N =", N)
  colnames(p_s) <- paste("m =", m)
  rownames(p_trho) <- paste("N =", N)
  colnames(p_trho) <- paste("m =", m)
  rownames(p_t3) <- paste("N =", N)
  colnames(p_t3) <- paste("m =", m)
  return(
    list(Snm = p_s, T_rho = p_trho, Tnm = p_t, T3 = p_t3)
  )
}


