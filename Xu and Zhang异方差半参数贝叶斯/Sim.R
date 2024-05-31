rm(list = ls())
setwd('D:\\Desktop\\Xu and Zhang异方差半参数贝叶斯')
source('Settings and Funcs.R')

n = 70
# Gibbs sampling
beta.res  <- matrix(0, nsim, 3)
gamma.res <- matrix(0, nsim, 3)
alpha.res <- matrix(0, nsim, K)
for(i in 1:nsim){
  cat('Sim = ', i, '\n', sep = "")
  para.process   <- Gibbs_sample(n, burn_in, mc_num)
  beta.process   <- para.process$beta.process
  gamma.process  <- para.process$gamma.process
  alpha.process  <- para.process$alpha.process
  beta_mc        <- tail(beta.process, mc_num)
  gamma_mc       <- tail(gamma.process, mc_num)
  alpha_mc       <- tail(alpha.process, mc_num)
  beta.res[i, ]  <- colMeans(beta_mc)
  gamma.res[i, ] <- colMeans(gamma_mc)
  alpha.res[i, ] <- colMeans(alpha_mc)
}
round(colMeans(beta.res), 3)
round(colMeans(gamma.res), 3)
round(apply(beta.res, 2, sd), 3)
round(apply(gamma.res, 2, sd), 3)

alpha.hat <- colMeans(alpha.res)
x.pre <- seq(0, 1, length.out = 100)
B.pre <- bs(x.pre, degree = 3, intercept = F)
y.pre <- B.pre %*% alpha.hat
plot(x.pre, 0.5*sin(2*pi*x.pre), type = 'l', ylim = c(-0.5, 0.5))
lines(x.pre, y.pre, col = "red", lwd = 2)

# para.process <- Gibbs_sample(n, burn_in, mc_num)
# u <- para.process$data$u
# plot(u, 0.5*sin(2*pi*u))
# plot(u, para.process$g.hat)




