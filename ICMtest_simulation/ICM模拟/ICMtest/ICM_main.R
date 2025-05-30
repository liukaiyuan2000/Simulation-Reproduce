rm(list =ls())
#设定新路径
setwd('D:/桌面/ICM模拟')
##  invoke libraries----
library(MASS)
source('ICM.R')


#######模拟
sample_size = c(rep(100, 4), 200, 400, 600)
dimensions <- c(2, 4, 6, 8, 12, 17, 20)
null_result <- 0
alternative_result <- 0
start <- Sys.time()
for(i in 1:length(dimensions)){
  null_result[i] <- ICMpower(sample_size[i], dimensions[i], 0)
  alternative_result[i] <- ICMpower(sample_size[i], dimensions[i], 0.25)
}
end <- Sys.time()
runningtime <- difftime(end, start, units = 'secs')
cat(runningtime, '\n')
null_result
alternative_result

write.csv(rbind(null_result, alternative_result), 
          'D:/桌面/ICM模拟结果/ICM-H11.csv')






