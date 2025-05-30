#' @description
#' Simulation codes for additive model. 
#' @author Meiling Chen, Kaiyuan Liu, Jiang Du 
#' Firstly, need to source functions code: source(<path>)
#' 
# source("D:/Desktop/Inv-BAM/funcs.R")
source("D:/Desktop/Inv-BAM/funcs.R")

n <- 100
################### Simulation Start ###################
gamma.res.c <- matrix(0, nsim, 2)
alpha.res.c <- map(1:p, \(x) matrix(0, nsim, K))
mu.res.c    <- rep(0, nsim)
gamma.res.r <- matrix(0, nsim, 2)
alpha.res.r <- map(1:p, \(x) matrix(0, nsim, K))
mu.res.r    <- rep(0, nsim)
rmseh.res   <- rep(0, nsim)
mae.res     <- rep(0, nsim)
ar.res      <- rep(0, nsim)
start.time  <- Sys.time()
for(i in 1:nsim){
  cat('Sim = ', i, '\n', sep = "")
  inv.data <- DGP(n, 'additive')
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
  mu.res.c[i]      <- mean(mu_mc.c)
  gamma.res.c[i, ] <- colMeans(gamma_mc.c)
  for(j in 1:p){
    alpha.res.c[[j]][i, ] <- map(alpha_mc.c, colMeans)[[j]]
  }
  
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
  mu.res.r[i]      <- mean(mu_mc.r)
  gamma.res.r[i, ] <- colMeans(gamma_mc.r)
  for(j in 1:p){
    alpha.res.r[[j]][i, ] <- map(alpha_mc.r, colMeans)[[j]]
  }
  
  c.hat <- mu.res.c[i] + rowSums(do.call(cbind, map2(alpha.res.c, inv.data$B.c, \(x, y) y %*% x[i, ])))
  r.hat <- mu.res.r[i] + rowSums(do.call(cbind, map2(alpha.res.r, inv.data$B.r, \(x, y) y %*% x[i, ])))
  l.hat <- c.hat - r.hat
  u.hat <- c.hat + r.hat
  l.tru <- inv.data$y.c - inv.data$y.r
  u.tru <- inv.data$y.c + inv.data$y.r
  
  rmseh.res[i] <- RMSEH(l.hat, l.tru, u.hat, u.tru)
  mae.res[i]   <- MAE(l.hat, l.tru, u.hat, u.tru)
  ar.res[i]    <- AR(l.hat, l.tru, u.hat, u.tru)
}
end.time <- Sys.time()
difftime(end.time, start.time)
######################### End ##########################
#' @export
#' 
mu.c
round(mean(mu.res.c), 3)
round(sd(mu.res.c), 3)
gamma.c
round(colMeans(gamma.res.c), 3)
round(apply(gamma.res.c, 2, sd), 3)
# mu.r
# round(mean(mu.res.r), 3)
# round(sd(mu.res.r), 3)
# gamma.r
# round(colMeans(gamma.res.r), 3)
# round(apply(gamma.res.r, 2, sd), 3)
round(mean(rmseh.res), 3)
round(sd(rmseh.res), 3)
round(mean(mae.res), 3)
round(sd(mae.res), 3)
round(mean(ar.res), 3)
round(sd(ar.res), 3)

{
  a1.hat <- map(alpha.res.c, colMeans)[[1]]
  x.pre <- seq(-3, 3, length.out = 100)
  B.pre <- scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
  y.pre <- B.pre %*% a1.hat
  plot(
    x.pre, x.pre, type = 'l', main = 'f1.c', 
    ylim = range(c(x.pre, y.pre))
  )
  lines(x.pre, y.pre, col = "red", lwd = 2)
  
  a2.hat <- map(alpha.res.c, colMeans)[[2]]
  x.pre <- seq(0, 1, length.out = 100)
  B.pre <-  scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
  y.pre <- B.pre %*% a2.hat
  plot(
    x.pre, sin(2*pi*x.pre), type = 'l', main = 'f2.c', 
    ylim = range(c(sin(2*pi*x.pre), y.pre))
  )
  lines(x.pre, y.pre, col = "red", lwd = 2)
  
  # a1.hat <- map(alpha.res.r, colMeans)[[1]]
  # x.pre <- seq(0, 5, length.out = 100)
  # B.pre <- scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
  # y.pre <- B.pre %*% a1.hat
  # plot(
  #   x.pre, x.pre, type = 'l', main = 'f1.r', 
  #   ylim = range(c(x.pre, y.pre))
  # )
  # lines(x.pre, y.pre, col = "red", lwd = 2)
  # 
  # a2.hat <- map(alpha.res.r, colMeans)[[2]]
  # x.pre <- seq(0, 1, length.out = 100)
  # B.pre <-  scale(bs(x.pre, degree = 3, knots = mean(x.pre)), scale = F)
  # y.pre <- B.pre %*% a2.hat
  # plot(
  #   x.pre, sin(2*pi*x.pre), type = 'l', main = 'f2.r', 
  #   ylim = range(c(sin(2*pi*x.pre), y.pre))
  # )
  # lines(x.pre, y.pre, col = "red", lwd = 2)
}





