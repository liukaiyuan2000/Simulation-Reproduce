#' @description
#' Simulation codes for additive model. \Parallel!!
#' @author Meiling Chen, Kaiyuan Liu, Jiang Du
#' Firstly, need to source functions code: source(<path>)
#'
# source("D:/Desktop/Inv-BAM/funcs.R")

source("D:/Desktop/Inv-BAM/funcs.R")

n = 100
no_cores = 6
nsim = 3*no_cores
case.choose = 'additive'
################### Simulation Start ###################
#' @example Gibbs sampling, \Parallel!!
#'
cl <- makeSOCKcluster(no_cores)
registerDoSNOW(cl)
start.time <- Sys.time()
RES <- foreach(
  i = 1:nsim, .options.snow = opts,
  .packages = c("splines", "MASS", "purrr"), .combine = bind_fun
) %dopar% {
  inv.data <- DGP(n, case.choose)

  para.process.c   <- Gibbs_sample(
    x = inv.data$x.c, z = inv.data$z.c, y = inv.data$y.c,
    B = inv.data$B.c, burn_in = burn_in, mc_num = mc_num
  )
  gamma.process.c  <- para.process.c$gamma.process
  alpha.process.c  <- para.process.c$alpha.process
  mu.process.c     <- para.process.c$mu.process
  gamma_mc.c       <- tail(gamma.process.c, mc_num)
  alpha_mc.c       <- map(alpha.process.c, \(x) tail(x, mc_num))
  mu_mc.c          <- tail(mu.process.c, mc_num)

  para.process.r   <- Gibbs_sample(
    x = inv.data$x.r, z = inv.data$z.r, y = inv.data$y.r,
    B = inv.data$B.r, burn_in = burn_in, mc_num = mc_num
  )
  gamma.process.r  <- para.process.r$gamma.process
  alpha.process.r  <- para.process.r$alpha.process
  mu.process.r     <- para.process.r$mu.process
  gamma_mc.r       <- tail(gamma.process.r, mc_num)
  alpha_mc.r       <- map(alpha.process.r, \(x) tail(x, mc_num))
  mu_mc.r          <- tail(mu.process.r, mc_num)

  mu.res.c = mean(mu_mc.c)
  gamma.res.c = colMeans(gamma_mc.c)
  alpha.res.c = map(alpha_mc.c, colMeans)
  mu.res.r = mean(mu_mc.r)
  gamma.res.r = colMeans(gamma_mc.r)
  alpha.res.r = map(alpha_mc.r, colMeans)
  c.hat <- mu.res.c + rowSums(do.call(cbind, map2(alpha.res.c, inv.data$B.c, \(x, y) y %*% x)))
  r.hat <- mu.res.r + rowSums(do.call(cbind, map2(alpha.res.r, inv.data$B.r, \(x, y) y %*% x)))
  l.hat <- c.hat - r.hat
  u.hat <- c.hat + r.hat
  l.tru <- inv.data$y.c - inv.data$y.r
  u.tru <- inv.data$y.c + inv.data$y.r

  rmseh.res <- RMSEH(l.hat, l.tru, u.hat, u.tru)
  mae.res   <- MAE(l.hat, l.tru, u.hat, u.tru)
  ar.res    <- AR(l.hat, l.tru, u.hat, u.tru)

  return(
    list(
      mu.res.c = mu.res.c, gamma.res.c = gamma.res.c,
      alpha.res.c = alpha.res.c,
      mu.res.r = mu.res.r, gamma.res.r = gamma.res.r,
      alpha.res.r = alpha.res.r,
      rmseh.res = rmseh.res, mae.res = mae.res, ar.res = ar.res
    )
  )
}

end.time <- Sys.time()
stopCluster(cl)
difftime(end.time, start.time)
######################### End ##########################
#' @export
#'
mu.c
round(mean(RES$mu.res.c), 3)
round(sd(RES$mu.res.c), 3)
gamma.c
round(colMeans(RES$gamma.res.c), 3)
round(apply(RES$gamma.res.c, 2, sd), 3)
round(mean(RES$rmseh.res), 3)
round(sd(RES$rmseh.res), 3)
round(mean(RES$mae.res), 3)
round(sd(RES$mae.res), 3)
round(mean(RES$ar.res), 3)
round(sd(RES$ar.res), 3)
# mu.r
# round(mean(RES$mu.res.r), 3)
# round(sd(RES$mu.res.r), 3)
# gamma.r
# round(colMeans(RES$gamma.res.r), 3)
# round(apply(RES$gamma.res.r, 2, sd), 3)
{
  if(case.choose == 'additive'){
    a1.hat <- map(RES$alpha.res.c, colMeans)[[1]]
    x.pre <- seq(-3, 3, length.out = 100)
    B.pre <- scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
    y.pre <- B.pre %*% a1.hat
    plot(
      x.pre, x.pre, type = 'l', main = 'f1.c',
      ylim = range(c(x.pre, y.pre))
    )
    lines(x.pre, y.pre, col = "red", lwd = 2)

    a2.hat <- map(RES$alpha.res.c, colMeans)[[2]]
    x.pre <- seq(0, 1, length.out = 100)
    B.pre <-  scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
    y.pre <- B.pre %*% a2.hat
    plot(
      x.pre, sin(2*pi*x.pre), type = 'l', main = 'f2.c',
      ylim = range(c(sin(2*pi*x.pre), y.pre))
    )
    lines(x.pre, y.pre, col = "red", lwd = 2)

    # a1.hat <- map(RES$alpha.res.r, colMeans)[[1]]
    # x.pre <- seq(0, 1, length.out = 100)
    # B.pre <- scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
    # y.pre <- B.pre %*% a1.hat
    # plot(
    #   x.pre, x.pre, type = 'l', main = 'f1.r',
    #   ylim = range(c(x.pre, y.pre))
    # )
    # lines(x.pre, y.pre, col = "red", lwd = 2)
    #
    # a2.hat <- map(RES$alpha.res.r, colMeans)[[2]]
    # x.pre <- seq(0, 1, length.out = 100)
    # B.pre <-  scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
    # y.pre <- B.pre %*% a2.hat
    # plot(
    #   x.pre, sin(2*pi*x.pre), type = 'l', main = 'f2.r',
    #   ylim = range(c(sin(2*pi*x.pre), y.pre))
    # )
    # lines(x.pre, y.pre, col = "red", lwd = 2)
  } else if(case.choose == 'linear'){
    a1.hat <- map(RES$alpha.res.c, colMeans)[[1]]
    x.pre <- seq(-3, 3, length.out = 100)
    B.pre <- scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
    y.pre <- B.pre %*% a1.hat
    plot(
      x.pre, 0.5*x.pre, type = 'l', main = 'f1.c',
      ylim = range(c(x.pre, y.pre))
    )
    lines(x.pre, y.pre, col = "red", lwd = 2)

    a2.hat <- map(RES$alpha.res.c, colMeans)[[2]]
    x.pre <- seq(-3, 3, length.out = 100)
    B.pre <-  scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
    y.pre <- B.pre %*% a2.hat
    plot(
      x.pre, x.pre, type = 'l', main = 'f2.c',
      ylim = range(c(sin(2*pi*x.pre), y.pre))
    )
    lines(x.pre, y.pre, col = "red", lwd = 2)
  }
}






