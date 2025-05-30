## Tnszsta: input(x, y) x，y -协变量和响应变量
##          output(Tnsz) -统计量
Tnszsta <- function(x, y){
  n = length(y)
  p = ncol(x)
  beta_n <- drop(solve(t(x)%*%x) %*% (t(x)%*%y))
  xnew <- drop(x %*% beta_n)
  
  x0 <- quantile(drop(xnew), 0.98)
  residu <- drop(y - xnew)
  
  ss <- outer(xnew, xnew, '<=')
  
  # sigma_n^2 --vector with dimensions n*1(sigma is a constant)
  sign2 <- mean(residu^2)
  
  # Rn1 --vector with dimensions n*1
  Rn1 <- colSums(ss * residu) / sqrt(n)
  
  # An --vector with dimension n*1
  alpha <- xnew / sign2
  An <- colMeans((alpha^2) * sign2 * t(ss))
  
  # Tnrn1 --vector with dimensions n*1
  sum1 <- rep(0, n)
  for(i in 1:n){
    temp1 <- matrix(ss[, i] * xnew / An, 1, n)
    temp2 <- matrix(residu * xnew, n, 1)
    sum1[i] <- temp1 %*% ss %*% temp2 / sign2 / sqrt(n) / n
  }
  Tnrn1 <- Rn1 - sum1
  
  ## test statistic(Wn) computation
  Wn <- mean((Tnrn1^2) * (xnew <= x0)) / 
    mean(xnew <= x0) / mean(xnew <= x0) / sign2
  return(Wn)
}





