source("D:/Desktop/学习/Bio2005/Funcs.R")
## H0: case = 1 VS H1: case = 2
cl <- makeSOCKcluster(no_cores)
registerDoSNOW(cl)
start.time <- Sys.time()
RES_sum <- Sim_fun(case = 2, m, n)

RES <- Res_fun(RES_sum)
end.time <- Sys.time()
stopCluster(cl)
difftime(end.time, start.time, units = 'secs')

RES