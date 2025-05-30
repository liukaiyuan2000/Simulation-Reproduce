## 计算矩阵的幂
"%^%" <- function(S, power) 
  with(eigen(S), vectors %*% (values ^ power * t(vectors)))

## CSE -- 计算B和qhat
CSE <- function(x, y){
  n <- nrow(x); p <- ncol(x);
  x_orig <- x
  x <- (x - matrix(rep(apply(x, 2, mean), n), n, p, byrow = T)) %*% (cov(x) %^% -0.5)
  ## P112-targetM(\hat{M})
  targetM <- array(0, c(p, p, n))
  for(i in 1:n){
    temp <- x * (y <= y[i])
    targetM[, , i] <- apply(temp, 2, mean) %*% t(apply(temp, 2, mean))
  }
  #求特征值和特征向量并按特征值从小到大排序
  fit <- eigen(as.matrix(apply(targetM, c(1, 2), mean)))
  lambda <- fit$values
  vectors <- fit$vectors
  directions <- (cov(x) %^% -0.5) %*% vectors
  directions <- directions[, 1:p]
  for(i in 1:p){
    directions[, i] <- directions[, i] / sqrt(sum(directions[, i]^2))
  }
  c = log(n) / n
  lambda1 = drop(t(lambda) ^ 2 + c)
  Q = c(lambda1[2:p], c) / lambda1[1:p]
  qhat = which(Q == min(Q))
  list(directions = directions, lambda = lambda, qhat = qhat, Q = Q)
}


 




