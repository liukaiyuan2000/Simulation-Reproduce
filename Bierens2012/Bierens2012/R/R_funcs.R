#' @title Data generate process
#'
#' @param Model choose model 1: linear model, 2: poisson model
#' @param case choose y case 1-4
#' @param n sample size, default = 200
#'
#' @return a list contains X and Y.
#' @export
#'
#' @examples data = DGP(Model = 1, case = 1)
#' x = data$X
#' y = data$Y
DGP <- function(Model, case, n = 200){
  x <- rnorm(n)
  switch(
    Model,
    {
      u <- switch(
        case,
        rnorm(n),
        rlogis(n),
        rt(n, 5),
        rnorm(n, 0, abs(x))
      )
      y <- 1 + x + u
    },
    {
      y <- switch(
        case,
        rpois(n, exp(x)),
        rnbinom(n, 1, 1 / (1+exp(-x))),
        rnbinom(n, 5, 1 / (1+exp(-x))),
        rnbinom(n, 10, 1 / (1+exp(-x)))
      )
    }
  )

  return(
    list(X = x, Y = y)
  )
}
#' @title Mapping Phi
#'
#' @param x a n*1 vector
#'
#' @return a one-to-one mapping Phi(x).
#' @export
#'
#' @examples x = rnorm(100)
#' phi.x = Phi(x)
Phi <- function(x){
  scale.x = (x - mean(x)) / sd(x)
  phi.x = atan(scale.x)
  return(phi.x)
}

#' @title Estimate model parameters
#'
#' @param x data X, n*1 vector
#' @param y data Y, n*1 vector
#' @param Model choose model 1: linear model, 2: poisson model
#'
#' @return the estimation of unknown parameters.
#' @export
#'
#' @examples data = DGP(Model = 1, case = 1)
#' x = data$X
#' y = data$Y
#' theta.hat = para.solve(x, y, Model = 1)
para.solve <- function(x, y, Model){
  switch(
    Model,
    {
      x = cbind(1, x)
      theta.hat <- solve(t(x) %*% x) %*% t(x) %*% y
    },
    {
      fit = glm(y~x, family = poisson(link = 'log'))
      theta.hat = summary(fit)$coef[, 1]
    }
  )
  return(theta.hat)
}

#' @title Sampling y.tilde
#'
#' @param x data X, n*1 vector
#' @param theta.hat the QMLE
#' @param Model choose model 1: linear model, 2: poisson model
#'
#' @return y.tilde, n*1 vector.
#' @export
#'
#' @examples data = DGP(Model = 1, case = 1)
#' x = data$X
#' y = data$Y
#' theta.hat = para.solve(x, y, Model = 1)
#' y.tilde = y.tilde.sample(x, theta.hat, Model = 1)
y.tilde.sample <- function(x, theta.hat, Model){
  n = length(x)
  y.tilde <- switch(
    Model,
    theta.hat[1] + theta.hat[2] * x + rnorm(n, 0, sd(x)),
    rpois(n, exp(theta.hat[1] + theta.hat[2] * x))
  )
  return(y.tilde)
}



#' @title Tnsc test statistic
#'
#' @param x data X, n*1 vector
#' @param y data Y, n*1 vector
#' @param y.tilde data Y.tilde, n*1 vector
#' @param c parameter c
#'
#' @return the value of Tnsc.
#' @export
#'
#' @examples data = DGP(Model = 1, case = 1)
#' x = data$X
#' y = data$Y
#' theta.hat = para.solve(x, y, Model = 1)
#' y.tilde = y.tilde.sample(x, theta.hat, Model = 1)
#' Tn.hat = Tnsc(x, y, y.tilde, c = 5)
#' @title Tnsc test statistic
#'
#' @param x data X, n*1 vector
#' @param y data Y, n*1 vector
#' @param y.tilde data Y.tilde, n*1 vector
#' @param c parameter c
#'
#' @return the value of Tnsc.
#' @export
#'
#' @examples data = DGP(Model = 1, case = 1)
#' x = data$X
#' y = data$Y
#' theta.hat = para.solve(x, y, Model = 1)
#' y.tilde = y.tilde.sample(x, theta.hat, Model = 1)
#' Tn.hat = Tnsc(x, y, y.tilde, c = 5)
Tnsc <- function(x, y, y.tilde, c){
  n = length(x)
  diffx <- outer(x, x, '-') * c
  diffy <- outer(y, y, '-') * c
  diffytilde <- outer(y.tilde, y.tilde, '-') * c
  diffyytilde <- outer(y, y.tilde, '-') * c
  diffytildey <- outer(y.tilde, y, '-') * c
  diffy.vec <- (y - y.tilde) * c
  I11 = sin(diffy) / diffy
  I12 = sin(diffytilde) / diffytilde
  I13 = sin(diffyytilde) / diffyytilde
  I14 = sin(diffytildey) / diffytildey
  I15 = sin(diffx) / diffx
  I11[is.na(I11)] = 1
  I12[is.na(I12)] = 1
  I13[is.na(I13)] = 1
  I14[is.na(I14)] = 1
  I15[is.na(I15)] = 1
  temp1 = (I11 + I12 - I13 - I14) * I15
  I1 = sum(temp1[upper.tri(temp1)])
  temp2 = sin(diffy.vec) / diffy.vec
  temp2[is.na(temp2)] = 1
  I2 = n - sum(temp2)
  return((I1 + I2) * (2 / n))
}




#' @title Bootstrap sample
#'
#' @param x data X, n*1 vector
#' @param theta.hat the QMLE
#' @param Model choose model 1: linear model, 2: poisson model
#' @param B the number of bootstrap
#'
#' @return y.star, a B*1 list.
#' @export
#'
#' @examples data = DGP(Model = 1, case = 1)
#' x = data$X
#' y = data$Y
#' theta.hat = para.solve(x, y, Model = 1)
#' y.tilde = y.tilde.sample(x, theta.hat, Model = 1)
#' Tn.hat = Tnsc(x, y, y.tilde, c = 5)
#' y.b = bootstrap.sample(x, theta.hat, Model = 1, B = 500)
bootstrap.sample <- function(x, theta.hat, Model, B){
  require(purrr)

  n = length(x)
  switch(
    Model,
    {
      sigma.hat = sd(x)
      u.star = map(1:B, \(x) rnorm(n, 0, sigma.hat))
      y.star = map(u.star, \(u) theta.hat[1] + theta.hat[2] * x + u)
    },
    {
      y.star = map(1:B, \(y) rpois(n, exp(theta.hat[1] + theta.hat[2] * x)))
    }
  )
  return(
    y.star = y.star
  )
}

#' Bootstrap test statistic sample
#'
#' @param x data X, n*1 vector
#' @param y.b bootstrap sample, B*1 list
#' @param c.vec the vector contains all parameter c
#' @param Model choose model 1: linear model, 2: poisson model
#' @param B the number of bootstrap
#' @param Phi.func the one-to-one mapping Phi
#'
#' @return a matrix contains: bootstrap test statistics in different c and max value.
#' @export
#'
#' @examples data = DGP(Model = 1, case = 1)
#' x = data$X
#' y = data$Y
#' theta.hat = para.solve(x, y, Model = 1)
#' y.tilde = y.tilde.sample(x, theta.hat, Model = 1)
#' Tn.hat = Tnsc(x, y, y.tilde, c = 5)
#' y.b = bootstrap.sample(x, theta.hat, Model = 1, B = 500)
#' Tn.b.hat = Tn.b(x, y.b, c.vec = 5, Model = 1, B = 500, Phi.func = NULL)
Tn.b <- function(x, y.b, c.vec, Model, B, Phi.func = NULL){
  require(purrr)
  theta.hat.b = map(y.b, \(y) para.solve(x, y, Model))
  y.tilde.b = map(theta.hat.b, \(theta) y.tilde.sample(x, theta, Model))

  if(is.null(Phi.func)){
    x = x
    y.b = y.b
  } else {
    x = Phi.func(x)
    y.b = map(y.b, Phi.func)
    y.tilde.b = map(y.tilde.b, Phi.func)
  }
  c.len <- length(c.vec)
  Tn.b = matrix(0, B, c.len)
  for(i in 1:c.len){
    Tn.b[, i] = map2_dbl(y.b, y.tilde.b, \(y, y.tilde) Tnsc_rcpp(x, y, y.tilde, c.vec[i]))
  }

  if(c.len > 1){
    Tn.b.max = apply(Tn.b, 1, max)
    return(cbind(Tn.b, Tn.b.max))
  } else{
    return(Tn.b)
  }

}







