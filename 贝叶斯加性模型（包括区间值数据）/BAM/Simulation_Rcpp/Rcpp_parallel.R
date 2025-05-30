#' @description
#' Simulation codes for additive model. \Parallel&\Rcpp!!
#' @author Meiling Chen, Kaiyuan Liu, Jiang Du
#' Firstly, need to source functions code: source(<path>)
#'
# source("D:/Desktop/陈美玲-加性异方差贝叶斯/additive_funs.R")

source("D:/Desktop/codes and literatures/陈美玲-加性异方差贝叶斯/Simulation_Rcpp/Rcpp_settings.R")

#' @param n sample size
#' @param no_cores number of cores, need to detect before setting (parallel::detectCores())
#'
n = 100
no_cores <- 5
#' @example Gibbs sampling, \Parallel!!
#'
cl <- makeSOCKcluster(no_cores)
registerDoSNOW(cl)
start.time <- Sys.time()
RES <- foreach(
  i = 1:nsim,.options.snow = opts, .packages = c("splines", "MASS", "purrr", "rBAM"), 
  .combine = bind_fun
) %dopar% {
  para.process   <- Gibbs_sample_rcpp(n, burn_in, mc_num)
  gamma.process  <- para.process$gamma.process
  alpha.process  <- para.process$alpha.process
  mu.process     <- para.process$mu.process
  gamma_mc       <- tail(gamma.process, mc_num)
  alpha_mc       <- map(alpha.process, \(x) tail(x, mc_num))
  mu_mc          <- tail(mu.process, mc_num)
  return(
    list(
      mu.res = mean(mu_mc), gamma.res = colMeans(gamma_mc),
      alpha.res = map(alpha_mc, colMeans)
    )
  )
}

end.time <- Sys.time()
stopCluster(cl)
difftime(end.time, start.time)
#' @export
#'
round(mean(RES$mu.res), 3)
round(colMeans(RES$gamma.res), 3)
round(sd(RES$mu.res), 3)
round(apply(RES$gamma.res, 2, sd), 3)

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






