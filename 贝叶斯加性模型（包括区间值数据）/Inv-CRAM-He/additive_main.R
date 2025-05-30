#' @description
#' Simulation codes for additive model. 
#' @author Meiling Chen, Kaiyuan Liu, Jiang Du 
#' Firstly, need to source functions code: source(<path>)
#' 
# source("D:/Desktop/陈美玲-加性异方差贝叶斯/additive_funs.R")

source("D:/Desktop/Inv-CRAM-He/additive_funs.R")

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
  y.true         <- para.process$y
  B              <- para.process$B
  gamma_mc       <- tail(gamma.process, mc_num)
  alpha_mc       <- map(alpha.process, \(x) tail(x, mc_num))
  mu_mc          <- tail(mu.process, mc_num)
  mu.res[i]      <- mean(mu_mc)
  gamma.res[i, ] <- colMeans(gamma_mc)
  m.hat <- 0
  for(j in 1:p){
    tempj <- map(alpha_mc, colMeans)[[j]]
    alpha.res[[j]][i, ] <- tempj
    m.hat <- m.hat + B[[j]] %*% tempj
  }
  y.hat <- m.hat + mu.res[i]
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
  a1.hat <- map(RES$alpha.res, colMeans)[[1]]
  x.pre <- seq(-1, 1, length.out = 100)
  B.pre <- scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
  y.pre <- B.pre %*% a1.hat
  plot(x.pre, x.pre, type = 'l', ylim = range(c(x.pre, y.pre)))
  lines(x.pre, y.pre, col = "red", lwd = 2)
}
{
  a2.hat <- map(RES$alpha.res, colMeans)[[2]]
  x.pre <- seq(0, 1, length.out = 100)
  B.pre <-  scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
  y.pre <- B.pre %*% a2.hat
  plot(x.pre, sin(2*pi*x.pre), type = 'l', ylim = range(c(sin(2*pi*x.pre), y.pre)))
  lines(x.pre, y.pre, col = "red", lwd = 2)
}
{
  a3.hat <- map(RES$alpha.res, colMeans)[[3]]
  x.pre <- seq(0, 1, length.out = 100)
  B.pre <-  scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
  y.pre <- B.pre %*% a3.hat
  plot(x.pre, cos(2*pi*x.pre), type = 'l', ylim = range(c(cos(2*pi*x.pre), y.pre)))
  lines(x.pre, y.pre, col = "red", lwd = 2)
}

su = 0
for (i in 1:length(B)) {
  s = B[[i]] %*% colMeans(alpha.res[[i]])
  su = s + su 
}
y.pre = su + mean(mu.res)

interval <- function(yc, yr){
  y_low = yc - 1/2*yr
  y_up = yc + 1/2*yr
  return(cbind(y_low, y_up))
}

y.true_r=y.true
y.pre_r=y.pre
y.true_c=y.true
y.pre_c=y.pre


y.interval_true = interval(y.true_c,y.true_r)
y.interval_pre = interval(y.pre_c,y.pre_r)
RMSE = sqrt((colMeans((y.interval_true - y.interval_pre)^2)))










