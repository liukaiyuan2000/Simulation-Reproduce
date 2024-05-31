rm(list = ls())
options(warn = -1)
library(doSNOW)
library(tcltk)
library(Bierens2012)
library(purrr)
nsim = 1000
no_cores <- 12
progress <- function(n){
  pb <- txtProgressBar(max = nsim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)

# For example: case = 1; Model = 1; alpha = 0.05; c = 5; B = 500
Bierens.test <- function(case, Model, alpha, c, B, nsim = 1000, Phi.func = NULL){
  res = rep(0, nsim)
  for(i in 1:nsim){
    cat(i, '\r')
    data = DGP(Model, case)
    x = data$X
    y = data$Y
    theta.hat = para.solve(x, y, Model)
    y.tilde = y.tilde.sample(x, theta.hat, Model)
    y.b = bootstrap.sample(x, theta.hat, Model, B)

    Tn.b.hat = Tn.b(x, y.b, c, Model, B, Phi.func)
    if(is.null(Phi.func)){
      Tn.hat = Tnsc_rcpp(x, y, y.tilde, c)
    } else {
      Tn.hat = Tnsc_rcpp(Phi.func(x), Phi.func(y), Phi.func(y.tilde), c)
    }

    res[i] = 1*(Tn.hat > sort(Tn.b.hat)[(1 - alpha)*B])
  }
  return(mean(res))
}

Bierens.test.parallel <- function(case, Model, alpha, c.vec, B = 500, nsim, Phi.func = NULL){
  c.len = length(c.vec)
  res <- foreach(
    sim = 1:nsim, .options.snow = opts, .packages = 'Bierens2012',
    .combine = 'rbind', .errorhandling = 'remove'
  )  %dopar% {
    data = DGP(Model, case)
    x = data$X
    y = data$Y
    theta.hat = para.solve(x, y, Model)
    y.tilde = y.tilde.sample(x, theta.hat, Model)
    y.b = bootstrap.sample(x, theta.hat, Model, B)
    Tn.b.hat = Tn.b(x, y.b, c.vec, Model, B, Phi.func)
    Tn.hat = rep(0, c.len)
    for(i in 1:c.len){
      if(is.null(Phi.func)){
        Tn.hat[i] = Tnsc_rcpp(x, y, y.tilde, c.vec[i])
      } else {
        Tn.hat[i] = Tnsc_rcpp(Phi.func(x), Phi.func(y), Phi.func(y.tilde), c.vec[i])
      }
    }
    Tn.hat = c(Tn.hat, max(Tn.hat))
    return(1*(Tn.hat > apply(Tn.b.hat, 2, sort)[(1 - alpha)*B, ]))
  }
  return(colMeans(res))
}

