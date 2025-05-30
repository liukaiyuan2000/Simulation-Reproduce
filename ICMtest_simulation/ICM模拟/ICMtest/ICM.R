## ICMpower: input(n, p, a) n-样本量；p-维度；a-参数(0或其他)
##           output(power) power-功效
# library(MASS)
ICMpower <- function(n, p, a){
  # s-模拟次数；c-bootstrap次数；result-结果储存
  s = 1000
  result = rep(0, s)
  c = 500
  mu = rep(0, p)
  ## case 1
  sigma = diag(rep(1, p)); c0 = rep(1, 2 * p);
  beta0 = rep(1, p) / sqrt(p)
  
  ## case 2
  # v = 1 / 2 ^ (0:(p - 1)); sigma = toeplitz(v);
  # p1 = floor(p / 2)
  # beta1 = c(rep(1, p1), rep(0, p - p1)) / sqrt(p1)
  # beta2 = c(rep(0, p - p1), rep(1, p1)) / sqrt(p1)
  for(l in 1:s){
    cat('模拟参数: ', 'p = ', p, '; n = ', n, '; a = ', a, '。第(', l, '/', 
        s, ')次循环', sep = "", '\r')
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
    residu <- drop(y - x %*% beta_n)
    
    x_new <- drop(x %*% beta_n)
    id1 <- rep(1:n, n)
    id2 <- rep(1:n, each = n)
    X <- x[id1, ] - x[id2, ]
    A <- matrix(apply(X, 1, crossprod), n, n)
    AA <- exp(-0.5 * A)
    # resi <- residu %*% t(residu)
    ## proposed test statistic
    # icm <- sum(AA * resi) / n
    icm <- drop(t(residu) %*% AA %*% residu / n)
    # 上一个公式的矩阵运算较慢，需要更改R中的dll文件
    
    # wild bootstrap
    simuicm <- rep(0, c)
    for(i in 1:c){
      prob <- (sqrt(5) - 1) / (2 * sqrt(5))
      z <- rbinom(n, 1, prob)
      u <- sqrt(5) * z + (1 - sqrt(5)) / 2
      ## estimation
      y_star <- x %*% beta_n + residu * u
      beta_star <- drop(solve(t(x)%*%x) %*% (t(x)%*%y_star))
      residu_star <- y_star - x %*% beta_star
      # resi_star <- residu_star %*% t(residu_star)
      # icm_star <- sum(AA * resi_star) / n
      icm_star <- drop(t(residu_star) %*% AA %*% residu_star / n)
      # 上一个公式的矩阵运算较慢，需要更改R中的dll文件
      simuicm[i] <- icm_star
    }
    simu <- sort(simuicm)
    critical <- simu[0.95 * c]
    result[l] <- 1 * (icm > critical)
  }
  cat('\nSimulation completed\n')
  power <- mean(result)
  return(power)
}















