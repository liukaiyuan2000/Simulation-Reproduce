rm(list =ls())
#设定新路径
setwd('D:/桌面/ICM模拟')
##  invoke libraries----
library(MASS)
source('Tnszpower.R')
source('BMfn.R')

Sim <-  5000
TN <- rep(0, Sim)
for(sim in 1:Sim){
  cat(sim,'\r')
  X <- BMfn(M = 500, T1 = 1)
  TN[sim] <- mean(X^2)
}
critical <- quantile(TN, .95) #1.647522
# critical <- 1.647522
#######模拟
sample_size <- c(rep(100, 4), 200, 400, 600)
dimensions <- c(2, 4, 6, 8, 12, 17, 20)
null_result <- 0
alternative_result <- 0
start <- Sys.time()
for(i in 1:length(dimensions)){
  null_result[i] <- Tnszpower(sample_size[i], dimensions[i], 0)
  alternative_result[i] <- Tnszpower(sample_size[i], dimensions[i], 0.1)
}
end <- Sys.time()
runningtime <- difftime(end, start, units = 'secs')
cat(runningtime, '\n')
null_result
alternative_result

write.csv(rbind(null_result, alternative_result), 'D:/桌面/ICM模拟结果/Tnsz-H14.csv')










