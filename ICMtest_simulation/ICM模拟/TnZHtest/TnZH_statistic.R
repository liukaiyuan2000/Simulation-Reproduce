## Tnzhsta: input(x, y, h) x，y -协变量和响应变量；h -带宽
##            output(TnZH) -统计量
TnZHsta <- function(x, y, h){
  beta_n <- drop(solve(t(x)%*%x) %*% (t(x)%*%y))
  xnew <- drop(x %*% beta_n)
  residu <- drop(y - xnew)
  
  id1 <- rep(1:n, n)
  id2 <- rep(1:n, each = n)
  X <- x[id1, ] - x[id2, ]
  A <- dnorm(sqrt(apply(X, 1, crossprod)) / h)
  AA <- matrix(A, n, n)
  diag(AA) <- 0
  temp1 <- t(residu) %*% AA %*% residu
  temp2 <- 2 * (t(residu)^2) %*% (AA^2) %*% (residu^2)
  ## test statistic computation
  TnZH <- drop(temp1 / sqrt(temp2))
  return(TnZH)
}





