rm(list =ls())
#设定新路径
setwd('D:/桌面/ICM模拟/杜老师的代码')
Rcpp::sourceCpp("SZtest.cpp")
Rcpp::sourceCpp("AICM.cpp")
Rcpp::sourceCpp("Escancianotest.cpp")
Rcpp::sourceCpp("ICM.cpp")
Rcpp::sourceCpp("Zhengtest.cpp")
# Rcpp::sourceCpp("ACM.cpp")

SIMpower <- function(n, p, a, b = 1000){
  # b-模拟次数；result-结果储存
  result = matrix(0, b, 5)
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
    tnzh <- Zhengtest(y, x, h)
    # aicm <- AICM(y, x)[2, ]
    aicm <- AICMsta(x, y)[2]
    tnsz <- SZtest(y, x)
    pcvm <- Escanciano(y, x)[2, ]
    icm <- ICM(y, x)[2, ]
    result[l, ] <- c(aicm, pcvm, tnsz, icm, tnzh)
  }
  cat('\nSimulation completed\n')
  result[, 1] <- 1 * (result[, 1] < 0.05)
  result[, 2] <- 1 * (result[, 2] < 0.05)
  result[, 3] <- 1 * (result[, 3] > 1.647522)
  result[, 4] <- 1 * (result[, 4] < 0.05)
  result[, 5] <- 1 * (result[, 5] > 1.96)
  power <- apply(result, 2, mean)
  return(power)
}

sample_size = c(rep(100, 4), 200, 400, 600)
dimensions <- c(2, 4, 6, 8, 12, 17, 20)
null_result <- matrix(0, 5, length(dimensions))
alternative_result <- matrix(0, 5, length(dimensions))
start <- Sys.time()
for(i in 1:length(dimensions)){
  null_result[, i] <- SIMpower(sample_size[i], dimensions[i], 0)
  alternative_result[, i] <- SIMpower(sample_size[i], dimensions[i], 0.25)
}
end <- Sys.time()
runningtime <- difftime(end, start, units = 'secs')
cat(runningtime, '\n')
null_result
alternative_result

write.csv(rbind(null_result, alternative_result), 'D:/桌面/ICM模拟结果/汇总.csv')

