Tntzpower <- function(n, p, a, b = 2000){
  # b-模拟次数；result-结果储存
  result = rep(0, b)
  ## case 1
  beta0 = rep(1, p) / sqrt(p)
  for(l in 1:b){
    cat('模拟参数: ', 'p = ', p, '; n = ', n, '; a = ', a, '。第(', l, '/', b, ')次循环', sep = "", '\r')
    x <- matrix(rnorm(n * p), n, p)
    ## H11
    y <- drop(x %*% beta0 + 0.25 * a * exp(x %*% beta0) + rnorm(n))
    ## test statistic computation
    tntz <- Tntzsta(x, y)
    result[l] <- tntz
  }
  cat('\nSimulation completed\n')
  result <- 1 * (result > 1.647522)
  power <- mean(result)
  return(power)
}
