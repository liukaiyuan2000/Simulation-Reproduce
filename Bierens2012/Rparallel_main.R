
source("D:/Desktop/Bierens2012/Rparallel_funcs.R")
Model = 2
c.vec = seq(5, 25, by = 5)
c.vec = 1:5
case.choose = 1:4
alpha.choose = c(0.01, 0.05, 0.1)

# tmp <- tempfile(fileext = ".out")
# Rprof(tmp)
cl <- makeSOCKcluster(no_cores)
registerDoSNOW(cl)
start.time <- Sys.time()
RES.sum <- c()
for(case in case.choose){
  for(alpha in alpha.choose){
    cat(
      '\nModel = ', switch(Model, 'Linear', 'Possion'),
      '; Case = ', case, '; Alpha = ', alpha,
      '\n', sep = ""
    )
    res <- Bierens.test.parallel(case, Model, alpha, c.vec, B = 500, nsim, Phi.func = NULL)
    RES.sum <- cbind(RES.sum, res)
  }
}

end.time <- Sys.time()
stopCluster(cl)
difftime(end.time, start.time, units = 'hours')
# Rprof(NULL)
# summaryRprof(tmp)

RES.sum

library(stargazer)
stargazer(RES.sum, summary = F, table.placement = "ht")



