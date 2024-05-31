library(splines)
library(MASS)
nsim = 20
burn_in = 100
mc_num = 200
beta = gamma = c(1, -0.5, 0.5)
# beta0 = gamma0 = rep(0, 3)
beta0 = gamma0 = c(1, -0.5, 0.5)
K = 6
alpha0 = rep(0, K)
b_beta = B_gamma = diag(3)
a_tau = b_tau = 1

# data generate process
DGP <- function(n){
  x <- matrix(runif(3*n, -1, 1), n, 3)
  z <- matrix(runif(3*n, -1, 1), n, 3)
  u <- runif(n)
  mu <- drop(x %*% beta + 0.5*sin(2*pi*u))
  sig2 <- drop(exp(z %*% gamma))
  eps <- rnorm(n)*sqrt(sig2)
  y <- mu + eps
  # y <- drop(x %*% beta) + eps
  B <- bs(u, degree = 3, knots = c(0.25, 0.5, 0.75))
  return(
    list(
      x = x, y = y, z = z, mu = mu, 
      u = u, B = B, sig2 = sig2, 
      eps = eps
    )
  )
}


p_gamma <- function(gamma, beta, alpha, x, y, z, B){
  Sig <- diag(drop(exp(z %*% gamma)))
  f <- det(Sig)^(-1/2)*
    exp(
      -1/2*(
        t(y - drop(x%*%beta) - drop(B%*%alpha)) %*% 
          solve(Sig) %*% (y - drop(x%*%beta) - drop(B%*%alpha)) + 
          t(gamma - gamma0) %*% solve(B_gamma) %*% (gamma - gamma0)
      )
    )
  return(drop(f))
}


Gibbs_sample <- function(n, burn_in, mc_num){
  data <- DGP(n)
  x <- data$x
  z <- data$z
  y <- data$y
  u <- data$u
  B <- data$B
  
  N = burn_in + mc_num
  beta.old  = beta0
  alpha.old = alpha0
  gamma.old = gamma0
  Sigma.old = diag(drop(exp(z %*% gamma.old)))
  
  beta.process  <- matrix(0, N, 3)
  alpha.process <- matrix(0, N, K)
  gamma.process <- matrix(0, N, 3)
  p_acc <- rep(0, N)
  
  for(s in 1:N){
    cat(s, '\r')
    
    # calculate tau2.new
    tau2.new <- 1/rgamma(1, K/2 + a_tau, 1/2*(t(alpha.old - alpha0) %*% (alpha.old - alpha0) + 2*b_tau))
    
    # calculate alpha.new
    b_alpha.star <- solve(1/tau2.new*diag(K) + t(B) %*% solve(Sigma.old) %*% B)
    alpha0.star  <- drop(b_alpha.star%*%(1/tau2.new*diag(K) %*% alpha0 + t(B) %*% solve(Sigma.old) %*% (y - x %*% beta.old)))
    alpha.new    <- mvrnorm(1, alpha0.star, b_alpha.star)
    # alpha.new <- alpha.old
    
    # calculate beta.new
    b_beta.star <- solve(solve(b_beta) + t(x) %*% solve(Sigma.old) %*% x)
    beta0.star  <- drop(b_beta.star %*% (solve(b_beta) %*% beta0 + t(x) %*% solve(Sigma.old) %*% (y - B %*% alpha.new)))
    beta.new    <- mvrnorm(1, beta0.star, b_beta.star)
    # beta.new = beta.old
    
    # calculate gamma.new, MH Algorithm
    Omega_gamma <- 1/2*t(z) %*% diag(drop((y - x%*%beta.new - B%*%alpha.new)^2/exp(z %*% gamma.old))) %*% z + solve(B_gamma)
    sig_gamma2  <- 1.5
    mc.new <- mvrnorm(1, gamma.old, sig_gamma2*solve(Omega_gamma))
    p_acc[s]  <- min(
      1, p_gamma(mc.new, beta.new, alpha.new, x, y, z, B) / 
        p_gamma(gamma.old, beta.new, alpha.new, x, y, z, B)
    )
    gamma.new <- gamma.old + (mc.new - gamma.old)*(runif(1) < p_acc[s])
    # gamma.new = gamma.old
    beta.process[s, ]  <- beta.old
    alpha.process[s, ] <- alpha.old
    gamma.process[s, ] <- gamma.old
    
    beta.old  = beta.new
    alpha.old = alpha.new
    gamma.old = gamma.new
    Sigma.old = diag(drop(exp(z %*% gamma.new)))
  }
  return(
    list(
      data = data, 
      gamma.process = gamma.process, 
      beta.process = beta.process, 
      alpha.process = alpha.process, 
      g.hat = drop(B %*% colMeans(tail(alpha.process, mc_num)))
    )
  )
}



