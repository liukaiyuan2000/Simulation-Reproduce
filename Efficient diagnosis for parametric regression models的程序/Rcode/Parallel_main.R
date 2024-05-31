source("D:/Desktop/Efficient diagnosis for parametric regression models的程序/Rcode/Parallel_funcs.R")
setwd("D:/Desktop/Efficient diagnosis for parametric regression models的程序/Rcode")
## H0: C = 0 VS H1: C != 0
setting.choose = 1:4
n.choose = c(100, 200, 300)
C.choose = cbind(0.2*(0:4), matrix(rep(0.1*(0:4), 3), 5, 3))

cl <- makeSOCKcluster(no_cores)
registerDoSNOW(cl)

RES.sum <- c()
start.time <- Sys.time()
for(setting in setting.choose){
  for(n in n.choose){
    for(C in C.choose[, setting]){
      cat('Simulation: setting = ', setting, ', n = ', n, ', C = ', C, '. \n', sep = "")
      RES.single <- Sim.fun(n, C, alpha, m, setting, nsim, rho)
      RES.sum <- rbind(RES.sum, RES.single)
      cat('\nDone~~\n')
    }
  }
}
end.time <- Sys.time()
stopCluster(cl)
difftime(end.time, start.time, units = 'mins')

RES.sum

library(stargazer)
stargazer(RES.sum, summary = F, table.placement = "H")

