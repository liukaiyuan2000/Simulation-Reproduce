rm(list=ls())
library(doSNOW)
library(tcltk)

Sim <- 2000
Sim.cv <- 3000
alpha <- 0.05
no_cores <- 12
progress.cv <- function(n){
  pb <- txtProgressBar(max = Sim.cv, style = 3, char = "->")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts.cv <- list(progress = progress.cv)

progress.sim <- function(n){
  pb <- txtProgressBar(max = Sim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts.sim <- list(progress = progress.sim)

critical.value <- function(n, p, Sim){
  
  ECP <- foreach(
    sim = 1:Sim, .options.snow = opts.cv, .packages = 'rRAT', .combine = 'rbind'
  )  %dopar% {
    X <- matrix(runif(n*p),n,p)
    Tn <- allmain(X)
    Pnp <- rep(0, 6)
    #sum-sum
    Pnp[1] <- 2*(1-pnorm(abs(Tn[1])))
    Pnp[2] <- 2*(1-pnorm(abs(Tn[2])))
    Pnp[3] <- 2*(1-pnorm(abs(Tn[3])))
    Pnp[4] <- 1-exp(-exp(-(Tn[4]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
    Pnp[5] <- 1-exp(-exp(-(Tn[5]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
    Pnp[6] <- 1-exp(-exp(-(Tn[6]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
    Cnp = min(Pnp)
    return(
      c(Tn[1:3], Tn[4:6]-4*log(p)+log(log(p)), Cnp)
    )
  }
  return(ECP)
}


Sim_func <- function(cvs, n, p, example.choose, case.choose, Sim){
  ECP <- foreach(
    sim = 1:Sim, .options.snow = opts.sim, .packages = c('rRAT', 'MASS'), .combine = 'rbind'
  )  %dopar% {
    X <- DGP(n, p, example.choose, case.choose)
    Pnp = rep(0, 6)
    Tn <- drop(allmain(X))

    Pnp[1] <- 2*(1-pnorm(abs(Tn[1])))
    Pnp[2] <- 2*(1-pnorm(abs(Tn[2])))
    Pnp[3] <- 2*(1-pnorm(abs(Tn[3])))
    Pnp[4] <- 1-exp(-exp(-(Tn[4]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
    Pnp[5] <- 1-exp(-exp(-(Tn[5]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
    Pnp[6] <- 1-exp(-exp(-(Tn[6]-4*log(p)+log(log(p)))/2)/(sqrt(8*pi)))
    Cnp = min(Pnp)
    return(
      c(
        1*((Tn[1]<cvs[1])||(Tn[1]>cvs[2])), 
        1*((Tn[2]<cvs[3])||(Tn[2]>cvs[4])), 
        1*((Tn[3]<cvs[5])||(Tn[3]>cvs[6])), 
        1*((Tn[4]-4*log(p)+log(log(p)))>cvs[7]), 
        1*((Tn[5]-4*log(p)+log(log(p)))>cvs[8]), 
        1*((Tn[6]-4*log(p)+log(log(p)))>cvs[9]), 
        1*(Cnp<cvs[10])
      )
    )
  }
  return(ECP)
}

