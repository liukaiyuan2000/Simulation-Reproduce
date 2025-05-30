source("D:/Desktop/张瑞云师姐/Funcs.R")
cl <- makeSOCKcluster(no_cores)
registerDoSNOW(cl)
start.time <- Sys.time()
RES <- Sim_ks(R, m, rho, case.DGP = 1, case.var = 1, alpha = 0.05, B = 500)

end.time <- Sys.time()
stopCluster(cl)
difftime(end.time, start.time, units = 'mins')

mean(RES)
