a <- seq(0, 1, by = 0.2)
n <- c(50, 100)
p = 8
result <- matrix(0, length(a), 2)
start <- Sys.time()
for(i in 1:length(a)){
  result[i, 1] <- Tntzpower(n[1], p, a[i])
  result[i, 2] <- Tntzpower(n[2], p, a[i])
}
end <- Sys.time()
runningtime <- difftime(end, start, units = 'secs')
cat(runningtime, '\n')
result


write.csv(rbind(null_result, alternative_result), 'D:/桌面/Tan2018test-H21.csv')


