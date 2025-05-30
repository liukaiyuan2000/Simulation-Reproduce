## Tntzsta: input(x, y) x，y -协变量和响应变量
##          output(Wn) -统计量
source("D:/桌面/ICM模拟/AICMtest/AICM_CSE().R")
Tntzsta <- function(x, y){
  n = length(y)
  p = ncol(x)
  beta_n <- drop(solve(t(x)%*%x) %*% (t(x)%*%y))
  directions <- CSE(x, y)$directions
  qhat <- CSE(x, y)$qhat
  # vec when qhat = 1;mat when qhat > 1
  B_n <- directions[, 1:qhat]
  alpha_hat <- rep(1, qhat) / sqrt(qhat)
  # n*1 vec
  xnew <- drop(x %*% B_n %*% alpha_hat)
  xnew1 <- drop(x %*% beta_n)
  
  x0 <- quantile(drop(xnew1), 0.99)
  residu <- drop(y - xnew)
  residu1 <- drop(y - xnew1)
  
  
  ss <- outer(xnew, xnew, '<=')

  # sigma_n^2 --vector with dimensions n*1(sigma is a constant)
  sign2 <- mean(residu1^2)

  # Vn --vector with dimensions n*1
  Vn <- colSums(ss * residu1) / sqrt(n)
  
  # An --vector with dimension n*1
  alpha <- xnew / sign2
  An <- colMeans((alpha^2) * sign2 * t(ss))
  
  # Tnrn1 --vector with dimensions n*1
  sum1 <- rep(0, n)
  for(i in 1:n){
    temp1 <- matrix(ss[, i] * xnew / An, 1, n)
    temp2 <- matrix(residu1 * xnew, n, 1)
    sum1[i] <- temp1 %*% ss %*% temp2 / sign2 / sqrt(n) / n
  }
  Tnrn1 <- Vn - sum1
  
  ## test statistic(Wn) computation
  Wn <- mean((Tnrn1^2) * (xnew1 <= x0)) / 
    mean(xnew1 <= x0) / mean(xnew1 <= x0) / sign2
  return(Wn)
}





