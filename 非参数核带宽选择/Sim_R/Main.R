rm(list = ls())
source("D:/Desktop/非参数核带宽选择/Funcs.R")
## Bandwidth Selection in Nonparametric Kernel Testing - JASA
## Only Proposed Method
n = 250
alpha = 0.95
Exp.case = 2
if(Exp.case == 1){
  Sim1.case = 3
  Sim2.case = NULL
  Sim2_e.case = NULL
  Sim3.case = NULL
  Sim4.case = NULL
  B <- 250
  kern.case = 1
} else if(Exp.case == 2){
  Sim1.case = NULL
  Sim2.case = 3
  Sim2_e.case = 3
  Sim3.case = NULL
  Sim4.case = NULL
  B <- switch(Sim2.case, 99, 250, 250)
  kern.case = 2
} else if(Exp.case == 3){
  Sim1.case = NULL
  Sim2.case = NULL
  Sim2_e.case = NULL
  Sim3.case = 2
  Sim4.case = NULL
  B <- 200
  kern.case = 1
} else{
  Sim1.case = NULL
  Sim2.case = NULL
  Sim2_e.case = NULL
  Sim3.case = NULL
  Sim4.case = 1
  B <- 200
  kern.case = 1
}
int_K_square <- ests_k(kern.case)$int_K_square
Kern_conv3_zero <- ests_k(kern.case)$Kern_conv3_zero

##### Start ----
registerDoSNOW(cl)
start.time <- Sys.time()
res <- foreach(
  sim = 1:nsim, .options.snow = opts, .packages = c('locfit', 'doSNOW', 'evd'), .combine = 'c'
)  %dopar% {
  data <- DGP(
    n, Example.case = Exp.case, Sim1.case, Sim2.case, Sim2_e.case, Sim3.case, Sim4.case
  )
  y <- data$y
  X <- data$X
  ests_x.store <- ests_x(X, Example.case = Exp.case)
  ests_yx.store <- ests_yx(y, X, Example.case = Exp.case)
  T_n.star.hat.restore <- rep(0, B)
  for(b in 1:B){
    e.star <- rnorm(n)
    y.star <- ests_yx.store$m.hat + sqrt(ests_yx.store$mu_2.hat)*e.star
    est.star <- ests_yx(y.star, X, Example.case = Exp.case)
    T_n.star.hat.restore[b] = T_n(y.star, X, ests_x.store, est.star)
  }
  l_alpha.star <- quantile(T_n.star.hat.restore, alpha)
  
  return(
    1*(T_n(y, X, ests_x.store, ests_yx.store) > l_alpha.star)
  )
}

end.time <- Sys.time()
stopCluster(cl)
difftime(end.time, start.time, units = 'secs')
##### End #####

mean(res)
