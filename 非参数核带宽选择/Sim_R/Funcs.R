##### invoke packages and some essential parameters ----
library(tcltk)
library(doSNOW)
library(locfit)
library(evd)

no_cores <- 16
nsim = 1000
progress <- function(n){
  pb <- txtProgressBar(max = nsim, style = 3, char = "*")
  f <- setTxtProgressBar(pb, n)
  return(f)
}
opts <- list(progress = progress)
cl <- makeSOCKcluster(no_cores)
#####
#' @description
#' Function of kernel function
#' 
#' @param x predictive variable
#' @param kern.case choose kernel function, 1-Gaussian, 2-Quartic
K <- function(x, kern.case){
  res <- switch(
    kern.case, 
    dnorm(x), 
    15/16*(1 - x^2)^2*(abs(x) <= 1)
  )
  return(res)
}
#' @description
#' Function of estimates about kernel function 
#' 
#' @param kern.case same as K()
ests_k <- function(kern.case){
  k = kern.case
  K3 <- function(x, y, z) {
    K(x, k)*K(y - x, k)*K(z - y, k)*K(-z, k)
  }
  K_square <- function(x){
    K(x, k)^2
  }
  int_K_square <- integrate(K_square, lower = -Inf, upper = Inf)$value
  Kern_conv3_zero <- integrate(function(z, y, x) {
    sapply(z, function(z) {
      integrate(function(y) {
        sapply(y, function(y) {
          integrate(function(x) {
            K3(x, y, z)
          }, -Inf, Inf)$value
        })
      }, -Inf, Inf)$value
    })
  }, -Inf, Inf)$value
  return(list(
    int_K_square = int_K_square, 
    Kern_conv3_zero = Kern_conv3_zero
  ))
}


#' @description 
#' Function of Data Generate Process
#' 
#' @param n sample size
#' @param Example.case choose one of the simulation examples 1 or 2
#' @param Sim1.case choose one of the cases 1-H0, 2-H1(j=1), and 3-H1(j=2) for Example 1
#' @param Sim2.case choose one of the cases 1-H0, 2-H1(tau=1), and 3-H1(tau=.25) for Example 2
#' @param Sim2_e.case choose one of the error cases 1-(a), 2-(b) and 3-(c) for Example 2
#' @param Sim3.case choose one of the cases 1-H0, 2-H1 for heteroskedasticity
#' @param Sim4.case choose one of the cases 1-H0, 2-H1 for e depend on X
DGP <- function(
    n, Example.case, Sim1.case = NULL, Sim2.case = NULL, Sim2_e.case = NULL, 
    Sim3.case = NULL, Sim4.case = NULL
  ){
  res <- switch(
    Example.case, 
    {
      e <- rnorm(n)
      x1 <- rep(0, n)
      x2 <- rep(0, n)
      for(i in 2:n){
        x1[i] <- x1[i - 1]*0.5 + rnorm(1)
        x2[i] <- x2[i - 1]*0.5 + rnorm(1)
      }
      cn <- switch(
        Sim1.case, 
        0, 
        n^(-0.5)*sqrt(log(log(n))), 
        n^(-7/18)
      )
      y <- x1 + x2 + e + cn*(x1^2 + x2^2)
      X <- cbind(x1, x2)
    }, 
    {
      L = qnorm(0.05, 0, 5)
      U = qnorm(0.95, 0, 5)
      q_simulation = runif(n)
      q_temp_simulation = pnorm(L, 0, 5) + q_simulation*(pnorm(U, 0, 5) - pnorm(L, 0, 5))
      x = qnorm(q_temp_simulation, 0, 5)
      e <- switch(
        Sim2_e.case, 
        rnorm(n, 0, 2), 
        c(rnorm(0.9*n, 0, 5), rnorm(0.1*n, 0, 1.4)), 
        rgev(n, -0.9000799, 1.55939)
      )
      cn <- switch(
        Sim2.case, 
        0, 
        5, 
        20
      )
      Delta <- switch(
        Sim2.case, 
        0, 
        dnorm(x), 
        dnorm(4*x)
      )
      y <- 1 + x + e + cn*Delta
      X <- matrix(x, n, 1)
    }, 
    {
      u <- runif(n)
      e <- rnorm(n, 0, sqrt(exp(2*u)))
      x <- rnorm(n)
      cn <- switch(
        Sim3.case, 
        0, 
        5
      )
      y <- x + e + cn*x^2
      X <- matrix(x, n, 1)
    }, 
    {
      x <- rnorm(n)
      e <- exp(0.5*x)
      cn <- switch(
        Sim4.case, 
        0, 
        5
      )
      y <- x + e + cn*x^2
      X <- matrix(x, n, 1)
    }
  )
  return(
    list(
      X = X, 
      y = y
    )
  )
}

#' @description
#' Function of generate some unknown parameter estimates only based on predictors
#' 
#' @param X design matrix
#' @param Example.case same as DGP()
ests_x <- function(X, Example.case){
  k <- switch(Example.case, 1, 2, 1, 1)
  d <- ncol(X)
  D_sum <- apply(X, 2, \(x) as.matrix(dist(x)))

  ## banwidth bcv
  bcv <- apply(X, 2, \(x) kdeb(x, meth = 'LSCV')/2.5)
  
  ## sum(K_h(Xi-Xj))
  kern.restore <- rep(1, n^2)
  for(j in 1:d){
    kern.restore <- kern.restore*K(D_sum[, j]/bcv[j], k)
  }
  kern_sum_mat <- matrix(kern.restore, n, n)
  ## pi.hat, v_2.hat in Eqs(11)
  pi.hat <- colMeans(kern_sum_mat) / prod(bcv)
  v_2.hat <- mean(pi.hat^2)
  return(
    list(
      d = d, 
      kern_sum_mat = kern_sum_mat, 
      pi.hat = pi.hat, 
      v_2.hat = v_2.hat
    )
  )
}

#' @description
#' Function of generate some unknown parameter estimates based on both predictors and response
#' 
#' @param y response variable
#' @param X design matrix
#' @param Example.case same as DGP()
ests_yx <- function(y, X, Example.case){
  switch(
    Example.case, 
    {
      theta.hat <- solve(t(X)%*%X)%*%t(X)%*%y
      m.hat <- drop(X %*% theta.hat)
    }, 
    {
      theta.hat <- solve(t(X)%*%X)%*%t(X)%*%(y - 1)
      m.hat <- drop(X %*% theta.hat) + 1
    }, 
    {
      theta.hat <- solve(t(X)%*%X)%*%t(X)%*%y
      m.hat <- drop(X %*% theta.hat)
    }, 
    {
      theta.hat <- solve(t(X)%*%X)%*%t(X)%*%y
      m.hat <- drop(X %*% theta.hat)
    }
  )
  e.hat <- y - m.hat
  mu_2.hat <- mean(e.hat^2)
  return(
    list(
      m.hat = m.hat, 
      e.hat = e.hat, 
      mu_2.hat = mu_2.hat
    )
  )
}
#' @description
#' Function of calculate h_ew
#' 
#' @param y response variable, n*1 vec. Contain y.star. 
#' @param X design matrix. In this case, contain x1 and x2, n*2 mat
#' @param data_ests_x ests_x(X, Example.case)
#' @param data_ests_yx ests_yx(y, X)
h_ew.hat <- function(y, X, data_ests_x, data_ests_yx){
  d <- data_ests_x$d
  kern_sum_mat <- data_ests_x$kern_sum_mat
  pi.hat <- data_ests_x$pi.hat
  v_2.hat <- data_ests_x$v_2.hat
  e.hat <- data_ests_yx$e.hat
  mu_2.hat <- data_ests_yx$mu_2.hat
  
  ## Delta_n in Eqs(36)
  delta_n <- colSums(kern_sum_mat*e.hat) / colSums(kern_sum_mat)
  ## \hat{C}_n^2 in Eqs(36)
  C_n_square <- mean(delta_n^2*pi.hat) / (mu_2.hat*sqrt(2*v_2.hat*int_K_square))
  ## t_n.hat in Eqs(36)
  t_n.hat <- n*C_n_square
  ## c_pi.hat
  c_pi.hat <- v_2.hat/mean(pi.hat)^3
  ## a_1.hat
  a_1.hat <- c_pi.hat*(sqrt(2)*Kern_conv3_zero)/(3*sqrt(int_K_square)^3)
  ## h_ew.hat in Eqs(36)
  h_ew.hat <- a_1.hat^(-1/(2*d))*t_n.hat^(-3/(2*d))
  return(h_ew.hat)
}


#' @description
#' Function of calculate test statistic T_n(h) in Eqs (11)
#' 
#' @param y response variable, n*1 vec. Contain y.star. 
#' @param X design matrix. In this case, contain x1 and x2, n*2 mat
#' @param data_ests_x ests_x(X, Example.case)
#' @param data_ests_yx ests_yx(y, X)
T_n <- function(y, X, data_ests_x, data_ests_yx){
  d <- data_ests_x$d
  kern_sum_mat <- data_ests_x$kern_sum_mat
  pi.hat <- data_ests_x$pi.hat
  v_2.hat <- data_ests_x$v_2.hat
  m.hat <- data_ests_yx$m.hat
  e.hat <- data_ests_yx$e.hat
  mu_2.hat <- data_ests_yx$mu_2.hat
  
  sigma_n_square <- 2*mu_2.hat^2*v_2.hat*int_K_square
  temp1 <- drop(t(e.hat)%*%kern_sum_mat%*%e.hat - t(e.hat)%*%diag(diag(kern_sum_mat))%*%e.hat)
  temp2 <- n*sqrt(h_ew.hat(y, X, data_ests_x, data_ests_yx)^d*sigma_n_square)
  t_n.hat <- temp1/temp2
  return(t_n.hat)
}
















