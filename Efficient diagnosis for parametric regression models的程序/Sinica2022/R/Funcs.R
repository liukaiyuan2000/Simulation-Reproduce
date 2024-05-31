#' @title Data generate process
#'
#' @param n sample size
#' @param C the difference bewteen H0 and H1
#' @param setting 1-4 represent model (5.1-5.4)
#'
#' @return a list contains X.tilde, Y.tilde, Z, and U.
#' @export
#'
#' @examples data = DGP(n = 10, C = 0, setting = 1)
DGP <- function(n, C, setting) {
  eps <- 0.15 * rnorm(n)
  U <- runif(n)

  switch(setting,
         {
           p = 2
           beta = 2:3
           X <- matrix(runif(n * p, 1, 2), ncol = p)
           X.tilde <- cbind((1 + 0.3 * cos(2 * pi * U)), (1 + 0.2 * (U^2 - 1/3))) * X
           Z <- NULL
           Y <- X %*% beta + C * exp(0.5 * X[,2]) + eps  # real response: n*1
           Y.tilde <- (1 + 0.2 * cos(2 * pi * U)) * Y  # distorted response: n*1
         },
         {
           p = 2
           beta = 1:2
           p = length(beta)
           X <- matrix(runif(n * p, 1, 2), ncol = p)
           X.tilde <- cbind((1 + 0.3 * cos(2 * pi * U)), (1 + 0.2 * (U^2 - 1/3))) * X
           Z <- NULL
           Y <- beta[1] + X[,1] * ((1 + X[,2])^beta[2]) + C * exp(0.5 * X[,2]) + eps
           Y.tilde <- (1 + 0.2 * cos(2 * pi * U)) * Y
         },
         {
           p = 5
           beta = rep(1, 5)
           p = length(beta)
           X <- matrix(runif(n * p, 1, 2), ncol = p)
           X.tilde <- cbind(
             (1 + 0.3 * cos(2 * pi * U)), (1 + 0.2 * (U^2 - 1/3)),
             (U + 1/2), (1 + 0.2 * (U^2 - 1/3)), (U^2 + 2/3)
            ) * X
           Z <- NULL
           Y <- X %*% beta + 2 * C * exp(0.5 * X[,2]) + eps
           Y.tilde <- (1 + 0.2 * cos(2 * pi * U)) * Y
         },
         {
           p = 6
           q = 4
           beta1 = c(rep(1, 4), rep(-1, 2))
           beta2 = rep(-1, 4)
           X <- matrix(runif(n * p, 1, 2), ncol = p)
           X.tilde <- cbind(
             (1 + 0.3 * cos(2 * pi * U)), (1 + 0.3 * cos(2 * pi * U)),
             (1 + 0.3 * cos(2 * pi * U)), (1 + 0.2 * (U^2 - 1/3)),
             (1 + 0.2 * (U^2 - 1/3)), (1 + 0.2 * (U^2 - 1/3))
           ) * X
           Z <- matrix(runif(n * q, 1, 2), ncol = q)
           Y <- X %*% beta1 + Z %*% beta2 + 0.1 * C * exp(X[,1] + X[,2]) + eps
           Y.tilde <- (1 + 0.2 * cos(2 * pi * U)) * Y
         }
  )

  return(list(X.tilde = X.tilde, Y.tilde = drop(Y.tilde), Z = Z, U = U))
}

#' @title Function of Beta.hat
#' @description
#' Obtain the Least Square estimate
#'
#'
#' @param X distorted covariates
#' @param Z observed covariates
#' @param Y distorted response
#' @param setting 1-4 represent model (5.1-5.4)
#'
#' @return estimates of beta and residual
#' @export
#'
#' @examples beta_e.hat = beta.hat.func(X.tilde, Z, Y.tilde, 1)
beta.hat.func <- function(X, Z, Y, setting) {
  if (setting == 2) {
    g <- function(X, beta){
      beta[1] + X[,1] * (1 + X[,2])^beta[2]
    }
    fit <- nls(
      Y ~ g(X, beta),
      start = list(beta = rep(0, 2))
    )
    beta.hat = summary(fit)$coef[, 1]
    e.hat <- Y - beta.hat[1] - X[,1] * (1 + X[,2])^beta.hat[2]
  } else if (setting == 4) {
    V <- cbind(X, Z)
    beta.hat <- drop(solve(t(V) %*% V + 1e-6, t(V) %*% Y))
    e.hat <- Y - V %*% beta.hat
  } else {
    beta.hat <- drop(solve(t(X) %*% X + 1e-6, t(X) %*% Y))
    e.hat <- Y - X %*% beta.hat
  }
  return(list(beta.hat = beta.hat, e.hat = as.matrix(e.hat)))
}


#' @title Calibrating function
#'
#' @param x a n*1 vector
#' @param u a n*1 vector
#' @param kern kernel parameters
#'
#' @return a n*1 calibrated vector
#' @export
#'
#' @examples Y.hat = cali.func(Y.tilde, U, kerns)
cali.func <- function(x, u, kern) {
  n <- length(x)
  s0 <- colMeans(kern)
  s1 <- colMeans(u * kern)
  s2 <- colMeans(u^2 * kern)
  temp1 <- colMeans(t(s2 - s1 * t(u)) * kern * x)
  temp2 <- s0 * s2 - s1^2
  idx0 <- (temp2 == 0)
  temp1[idx0] <- 1
  temp2[idx0] <- 1
  Ex <- temp1 / temp2
  x.hat <- mean(x) * x / Ex
  return(x.hat)
}


#' @title Calculate Tcvm
#'
#' @param X X.hat
#' @param Z
#' @param e e.hat
#'
#' @return Tcvm test statistics(2.4)
#' @export
#'
#' @examples Tcvm(X.hat, Z, e.hat)
Tcvm <- function(X, Z, e){
  V <- cbind(X, Z)
  n <- nrow(V)
  J <- ncol(V)
  p <- ncol(X)
  C_p <- pi^(J/2 - 1) / gamma(J/2 + 1)
  A <- matrix(0, nrow = n, ncol = n)
  diag(A) = C_p * (n + 1) * pi
  temp <- list()
  for (i in 1:n) {
    temp[[i]] <- matrix(rep(X[i,], n), ncol = p, byrow = TRUE) - X
  }
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      temp1 <- rowSums(temp[[i]] * temp[[j]])
      temp2 <- sqrt(rowSums(temp[[i]]^2)) * sqrt(rowSums(temp[[j]]^2))
      temp1[temp2 == 0] <- 1
      temp2[temp2 == 0] <- 1
      temp3 = round(temp1 / temp2, 7)
      A[i, j] <- C_p * sum(abs(pi - acos(temp3)))
      A[j, i] <- A[i, j]
    }
  }

  Tn <- t(e) %*% A %*% e / n^2
  return(list(Tn = drop(Tn), A = A))
}

#' @title Calculate Tcvm.rand
#'
#' @param V cbind(X.hat, Z) a n*(p+q) matrix
#' @param theta
#' @param e e.hat
#'
#' @return Tcvm test statistics(2.5)
#' @export
#'
#' @examples Tcvm.rand(V, theta.normal, e.hat)
Tcvm.rand <- function(V, theta, e) {
  n   <- nrow(e)
  rho <- ncol(e)
  M   <- ncol(theta)
  Tn  <- matrix(0, nrow = M, ncol = rho)
  Xt  <- V %*% theta
  for (m in 1:M) {
    II <- outer(Xt[, m], Xt[, m], "<=")
    Tn[m, ] <- colSums((t(II) %*% e)^2) / n^2
  }
  Tn <- colMeans(Tn)
  return(Tn)
}

#' @title Calculate Tks.rand
#'
#' @param V cbind(X.hat, Z) a n*(p+q) matrix
#' @param theta
#' @param e e.hat
#'
#' @return Tks test statistics(3.1)
#' @export
#'
#' @examples Tks.rand(V, theta.normal, e.hat)
Tks.rand <- function(V, theta, e){
  n   <- nrow(e)
  rho <- ncol(e)
  M   <- ncol(theta)
  Xt  <- V %*% theta
  Bn  <- matrix(0, n*M, rho)
  for(i in 1:(n*M)){
    Bn[i, ] <- rowSums(abs(t(e) %*% (Xt <= Xt[i])))
  }
  Tn <- apply(Bn, 2, max) / (M*sqrt(n))
  return(Tn)
}

#' @title Bootstrap procedure
#'
#'
#' @param X.hat
#' @param Z
#' @param Y.hat
#' @param e.hat
#' @param setting
#' @param rho the number of Bootstrap replications
#'
#' @return the Bootstrap residuals, a n*rho matrix
#' @export
#'
#' @examples e.star = bootstrap.res(X.hat, Z, Y.hat, e.hat, 1, rho = 1000)
bootstrap.res <- function(X.hat, Z, Y.hat, e.hat, setting, rho) {
  e.hat <- drop(e.hat)
  n <- nrow(X.hat)
  e <- matrix(rnorm(n*rho), n, rho)
  Y.star <- Y.hat - e.hat + e.hat * e
  e.star <- apply(Y.star, 2, \(x) beta.hat.func(X.hat, Z, x, setting)$e.hat)
  return(e.star)
}






