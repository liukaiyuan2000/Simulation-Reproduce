#' @description
#' Simulation codes for additive model. 
#' @author Meiling Chen, Kaiyuan Liu, Jiang Du 
#' Firstly, need to source functions code: source(<path>)
#' 
# source("D:/Desktop/陈美玲-加性异方差贝叶斯/additive_funs.R")

source("D:/Desktop/codes and literatures/陈美玲-加性异方差贝叶斯/Simulation_R/additive_funs.R")

#' @param n sample size
#' @param gamma.res result for \gamma
#' @param mu.res result for \mu
#' @param alpha.res result for \alpha
#' 
n         <- 100
gamma.res <- matrix(0, nsim, 2)
alpha.res <- map(1:p, \(x) matrix(0, nsim, K))
mu.res    <- rep(0, nsim)
#' @example Gibbs sampling
#' 
start.time <- Sys.time()
for(i in 1:nsim){
  cat('Sim = ', i, '\n', sep = "")
  para.process   <- Gibbs_sample(n, burn_in, mc_num)
  gamma.process  <- para.process$gamma.process
  alpha.process  <- para.process$alpha.process
  mu.process     <- para.process$mu.process
  gamma_mc       <- tail(gamma.process, mc_num)
  alpha_mc       <- map(alpha.process, \(x) tail(x, mc_num))
  mu_mc          <- tail(mu.process, mc_num)
  mu.res[i]      <- mean(mu_mc)
  gamma.res[i, ] <- colMeans(gamma_mc)
  for(j in 1:p){
    alpha.res[[j]][i, ] <- map(alpha_mc, colMeans)[[j]]
  }
}
end.time <- Sys.time()
difftime(end.time, start.time)
#' @export
#' 
round(mean(mu.res), 3)
round(colMeans(gamma.res), 3)
round(sd(mu.res), 3)
round(apply(gamma.res, 2, sd), 3)

{
  a1.hat <- map(alpha.res, colMeans)[[1]]
  x.pre <- seq(-1, 1, length.out = 100)
  B.pre <- scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
  y.pre <- B.pre %*% a1.hat
  plot(x.pre, x.pre, type = 'l', ylim = range(c(x.pre, y.pre)))
  lines(x.pre, y.pre, col = "red", lwd = 2)
}
{
  a2.hat <- map(alpha.res, colMeans)[[2]]
  x.pre <- seq(0, 1, length.out = 100)
  B.pre <-  scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
  y.pre <- B.pre %*% a2.hat
  plot(x.pre, sin(2*pi*x.pre), type = 'l', ylim = range(c(sin(2*pi*x.pre), y.pre)))
  lines(x.pre, y.pre, col = "red", lwd = 2)
}
{
  a3.hat <- map(alpha.res, colMeans)[[3]]
  x.pre <- seq(0, 1, length.out = 100)
  B.pre <-  scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
  y.pre <- B.pre %*% a3.hat
  plot(x.pre, cos(2*pi*x.pre), type = 'l', ylim = range(c(cos(2*pi*x.pre), y.pre)))
  lines(x.pre, y.pre, col = "red", lwd = 2)
}





