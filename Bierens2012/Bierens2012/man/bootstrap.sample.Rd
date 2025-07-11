% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_funcs.R
\name{bootstrap.sample}
\alias{bootstrap.sample}
\title{Bootstrap sample}
\usage{
bootstrap.sample(x, theta.hat, Model, B)
}
\arguments{
\item{x}{data X, n*1 vector}

\item{theta.hat}{the QMLE}

\item{Model}{choose model 1: linear model, 2: poisson model}

\item{B}{the number of bootstrap}
}
\value{
y.star, a B*1 list.
}
\description{
Bootstrap sample
}
\examples{
data = DGP(Model = 1, case = 1)
x = data$X
y = data$Y
theta.hat = para.solve(x, y, Model = 1)
y.tilde = y.tilde.sample(x, theta.hat, Model = 1)
Tn.hat = Tnsc(x, y, y.tilde, c = 5)
y.b = bootstrap.sample(x, theta.hat, Model = 1, B = 500)
}
