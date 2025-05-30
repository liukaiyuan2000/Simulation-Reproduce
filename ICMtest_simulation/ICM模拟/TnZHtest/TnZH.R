## TnZHpower: input(n, p, a) n-样本量；p-维度；a-参数(0或其他)
##            output(power) power-功效
TnZHpower <- function(n, p, a){
  # b-模拟次数；result-结果储存
  b = 1000
  result = rep(0, b)
  mu = rep(0, p)
  ## case 1
  sigma = diag(p)
  beta0 = rep(1, p) / sqrt(p)
  
  ## case 2
  # v = 1 / 2 ^ (0:(p - 1)); sigma = toeplitz(v);
  # p1 = floor(p / 2)
  # beta1 = c(rep(1, p1), rep(0, p - p1)) / sqrt(p1)
  # beta2 = c(rep(0, p - p1), rep(1, p1)) / sqrt(p1)
  for(l in 1:b){
    cat('模拟参数: ', 'p = ', p, '; n = ', n, '; a = ', a, '。第(', l, 
        '/', b, ')次循环', sep = "", '\r')
    x <- mvrnorm(n, mu, sigma)
    ## H11
    y <- drop(x %*% beta0 + a * exp(-(x %*% beta0)^2) + rnorm(n))
    ## H12
    # y <- drop(x %*% beta0 + a * cos(0.6 * pi * x %*% beta0) + rnorm(n))
    ## H13
    # y <- drop(x %*% beta1 + a*(x %*% beta2)^2 + rnorm(n))
    ## H14
    # y = drop(x %*% beta1 + a * exp(x %*% beta2) + rnorm(n))
    
    beta_n <- drop(solve(t(x)%*%x) %*% (t(x)%*%y))
    xnew <- drop(x %*% beta_n)
    residu <- drop(y - xnew)
    #h <- n^(-0.4)
    h <- 1.5 * n^(-1 / (4 + p))
    
    id1 <- rep(1:n, n)
    id2 <- rep(1:n, each = n)
    X <- (x[id1, ] - x[id2, ])^2
    A <- dnorm(sqrt(apply(X, 1, sum)) / h)
    AA <- matrix(A, n, n)
    diag(AA) <- 0
    temp1 <- t(residu) %*% AA %*% residu
    temp2 <- 2 * (t(residu)^2) %*% (AA^2) %*% (residu^2)
    ## test statistic computation
    TnZH <- drop(temp1 / sqrt(temp2))
    ## critical 1.96
    result[l] <- 1 * (TnZH > 1.96)
  }
  cat('\nSimulation completed\n')
  power <- mean(result)
  return(power)
}





