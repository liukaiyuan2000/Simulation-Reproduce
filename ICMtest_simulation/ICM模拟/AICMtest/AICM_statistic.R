## AICMsta: input(x, y) n-样本量；p-维度；a-参数(0或其他)
##            output(aicm) -统计量
source("D:/桌面/ICM模拟/AICM_CSE().R")
AICMsta <- function(x, y){
  n = length(y)
  p = ncol(x)
  beta_n <- drop(solve(t(x)%*%x) %*% (t(x)%*%y))
  residu <- drop(y - x %*% beta_n)
  
  directions <- CSE(x, y)$directions
  qhat <- CSE(x, y)$qhat
  B_n <- directions[, 1:qhat]
  x_new1 <- x %*% B_n
  x_new2 <- drop(x %*% beta_n)
  id1 <- rep(1:n, n)
  id2 <- rep(1:n, each = n)
  temp1 <- as.matrix(x_new1[id1, ] - x_new1[id2, ], n * n, qhat)
  o1 <- matrix(apply(temp1, 1, crossprod), n, n)
  o2 <- outer(x_new2, x_new2, '-')^2
  A <- o1 + o2
  AA <- exp(-0.5 * A)
  resi <- residu %*% t(residu)
  ## proposed test statistic
  aicm <- sum(AA * resi) / n
  # aicm <- t(residu) %*% AA %*% residu / n
  # 上一个公式的矩阵运算较慢，需要更改R中的dll文件
  # wild bootstrap
  c = 500
  simuicm <- rep(0, c)
  prob <- (sqrt(5) - 1) / (2 * sqrt(5))
  for(i in 1:c){
    z <- rbinom(n, 1, prob)
    u <- sqrt(5) * z + (1 - sqrt(5)) / 2
    ## estimation
    y_star <- x %*% beta_n + residu * u
    beta_star <- drop(solve(t(x)%*%x) %*% (t(x)%*%y_star))
    residu_star <- y_star - x %*% beta_star
    resi_star <- residu_star %*% t(residu_star)
    aicm_star <- sum(AA * resi_star) / n
    # aicm_star <- drop(t(residu_star) %*% AA %*% residu_star / n)
    # 上一个公式的矩阵运算较慢，需要更改R中的dll文件
    simuicm[i] <- aicm_star
  }
  p_value <- mean(1 * (simuicm > aicm))
  return(c(aicm, p_value))
}















