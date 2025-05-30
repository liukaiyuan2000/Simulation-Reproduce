source("D:/Desktop/JASA-Qin/CL2011_func.R")
#####　Simulation for Model 1 with hard thresholding
p.choose = c(30, 100, 200)
n = 100
nsim = 500
model = 2

RES.mean = NULL
RES.sd = NULL

for(norm_type in c("operator", "l1", "frobenius")){
  for(p in p.choose){
    for(s_type in c("adaptive lasso", "hard")){
      for(method in c('universal', 'adaptive', 'adaptive2')){
        risk = numeric(nsim)
        for(i in 1:nsim){
          cat(s_type, norm_type, p, method, i, '\r', sep = "; ")
          mu = rep(0, p)
          Sigma0 = gen_sigma(p, model = model)
          X = mvrnorm(n, mu, Sigma0)
          Sigma.hat = cov(X)
          thre.cv = thre_cv_func(X, method = method, s_type = s_type)
          Sigma.cv = thresholding(Sigma.hat, thre.cv, s_type = s_type)
          risk[i] <- cal_risk(Sigma0, Sigma.cv, norm_type)
        }
        cat('\n    Done!\n')
        RES.mean = c(RES.mean, mean(risk))
        RES.sd = c(RES.sd, sd(risk) / sqrt(n))
      }
    }
  }
}

RES <- paste(
  sprintf('%.2f', RES.mean), "(", 
  sprintf('%.2f', RES.sd), ")", sep = ""
)
RES.dat = data.frame(matrix(RES, ncol = 6, byrow = T))
colnames(RES.dat) = c(
  'A-Universal', 'A-Adaptive', 'A-Adaptive2', 
  'H-Universal', 'H-Adaptive', 'H-Adaptive2'
)
rownames(RES.dat) = c(
  'O-30', 'O-100', 'O-200', 
  'L1-30', 'L1-100', 'L1-200', 
  'F-30', 'F-100', 'F-200'
)
RES.dat



