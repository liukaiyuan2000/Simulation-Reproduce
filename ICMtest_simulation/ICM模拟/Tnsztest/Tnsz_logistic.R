## Tnszpower: input(n) n-样本量
##            output(power) power-功效
Tnszpower_logistic <- function(n = 100, critical = 1.647522){
  # b-模拟次数；result-结果储存
  b = 2000
  result = rep(0, b)
  mu = c(0, 0, 0)
  sigma = diag(3)
  beta0 = c(-1, -1, 2)
  
  for(l in 1:b){
    cat('模拟参数: ', 'n = ', n, '。第(', l, '/', b, ')次循环', sep = "", '\r')
    x <- mvrnorm(n, mu, sigma)
    xx <- x %*% beta0
    # P <- drop(exp(xx) / (1 + exp(xx))) # null
    P <- pnorm(xx) # alternative
    y <- 1 * (runif(n) > P)

    beta_n <- glm(y ~ x + 0, family = "binomial")$coef
    xnew <- drop(x %*% beta_n)
    
    x0 <- quantile(drop(xnew), 0.98)
    residu <- drop(y - exp(xnew) / (1 + exp(xnew)))
    
    ss <- outer(xnew, xnew, '<=')
    
    # sigma_n^2 --vector with dimensions n*1(sigma is a constant)
    sign2 <- exp(xnew) / (1 + exp(xnew))^2
    
    # Rn1 --vector with dimensions n*1
    Rn1 <- colSums(ss * residu) / sqrt(n)
    
    # An --vector with dimension n*1
    An <- colMeans((xnew^2) * sign2 * t(ss))
    
    # Tnrn1 --vector with dimensions n*1
    sum1 <- rep(0, n)
    for(i in 1:n){
      temp1 <- matrix(ss[, i] * xnew * sign2 / An, 1, n)
      temp2 <- matrix(residu * xnew, n, 1)
      sum1[i] <- temp1 %*% ss %*% temp2 / sqrt(n) / n
    }
    Tnrn1 <- Rn1 - sum1
    
    ## test statistic(Wn) computation
    Wn <- mean((Tnrn1^2) * sign2 * (xnew <= x0)) / (mean(residu^2*(xnew<x0)))^2
    
    ## critical BMfn.R
    result[l] <- 1 * (Wn > critical)
  }
  power <- mean(result)
  return(power)
}





