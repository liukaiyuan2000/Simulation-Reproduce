rm(list =ls())
#设定新路径
setwd('D:/桌面/ICM模拟/杜老师的代码')
Rcpp::sourceCpp("Zhengtest.cpp")

TnZHpower <- function(n, p, a, b = 1000){
  # b-模拟次数；result-结果储存
  result = rep(b, 0)
  ## case 1
  beta0 = rep(1, p) / sqrt(p)
  for(l in 1:b){
    cat('模拟参数: ', 'p = ', p, '; n = ', n, '; a = ', a, '。第(', l, '/', b, ')次循环', sep = "", '\r')
    x <- matrix(rnorm(n * p), n, p)
    ## H11
    y <- drop(x %*% beta0 + a * exp(-(x %*% beta0)^2) + rnorm(n))
    # h <- n^(-0.4)
    h <- 1.5 * n^(-1 / (4 + p))
    ## test statistic computation
    # TnZH
    Tn <- Zhengtest(y, x, h)
    result[l] <- 1 * (Tn > 1.96)
  }
  cat('\nSimulation completed\n')
  power <- mean(result)
  return(power)
}

sample_size = c(rep(100, 4), 200, 400, 600)
dimensions <- c(2, 4, 6, 8, 12, 17, 20)
null_result <- 0
alternative_result <- 0
start <- Sys.time()
for(i in 1:length(dimensions)){
  null_result[i] <- TnZHpower(sample_size[i], dimensions[i], 0)
  alternative_result[i] <- TnZHpower(sample_size[i], dimensions[i], 0.25)
}
end <- Sys.time()
runningtime <- difftime(end, start, units = 'secs')
cat(runningtime, '\n')
null_result
alternative_result

write.csv(rbind(null_result, alternative_result), 'D:/桌面/ICM模拟结果/TnZH-H11.csv')








