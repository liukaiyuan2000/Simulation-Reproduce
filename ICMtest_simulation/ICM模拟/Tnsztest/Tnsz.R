## Tnszpower: input(n, p, a) n-样本量；p-维度；a-参数(0或其他)
##            output(power) power-功效
Tnszpower <- function(n, p, a, critical = 1.647522){
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
    
    ## critical BMfn.R
    result[l] <- 1 * (Wn > critical)
  }
  cat('\nSimulation completed\n')
  power <- mean(result)
  return(power)
}





