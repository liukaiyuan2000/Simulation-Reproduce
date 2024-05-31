source("D:/Desktop/CSTM2017/Funcs.R")
## H0: case = 1-3 VS H1: case = 4
cl <- makeSOCKcluster(no_cores)
registerDoSNOW(cl)
start.time <- Sys.time()
RES_sum <- Sim_fun(case = 1, m, N)

RES <- Res_fun(RES_sum)
end.time <- Sys.time()
stopCluster(cl)
difftime(end.time, start.time, units = 'secs')

RES