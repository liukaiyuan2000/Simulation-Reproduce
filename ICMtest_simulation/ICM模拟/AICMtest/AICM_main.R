rm(list =ls())
#设定新路径
setwd('D:/桌面/学习/ICMtest_simulation/ICM模拟/AICMtest')
##  invoke libraries----
library(MASS)
source('AICM_CSE().R')
source('AICM.R')


#######模拟
sample_size = c(rep(100, 2), 200)
dimensions <- c(2, 6, 12)
null_result <- 0
alternative_result <- 0
start <- Sys.time()
for(i in 1:length(dimensions)){
  null_result[i] <- AICMpower(sample_size[i], dimensions[i], 0)
  alternative_result[i] <- AICMpower(sample_size[i], dimensions[i], 0.25)
}
end <- Sys.time()
runningtime <- difftime(end, start, units = 'secs')
cat(runningtime, '\n')
null_result
alternative_result

write.csv(rbind(null_result, alternative_result), 
          'D:/桌面/ICM模拟结果/AICM-H12.csv')






