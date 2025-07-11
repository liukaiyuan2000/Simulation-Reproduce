% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_funcs.R
\name{Tnsc}
\alias{Tnsc}
\title{Tnsc test statistic}
\usage{
Tnsc(x, y, y.tilde, c)
}
\arguments{
\item{x}{data X, n*1 vector}

\item{y}{data Y, n*1 vector}

\item{y.tilde}{data Y.tilde, n*1 vector}

\item{c}{parameter c}
}
\value{
the value of Tnsc.

the value of Tnsc.
}
\description{
Tnsc test statistic

Tnsc test statistic
}
\examples{
data = DGP(Model = 1, case = 1)
x = data$X
y = data$Y
theta.hat = para.solve(x, y, Model = 1)
y.tilde = y.tilde.sample(x, theta.hat, Model = 1)
Tn.hat = Tnsc(x, y, y.tilde, c = 5)
data = DGP(Model = 1, case = 1)
x = data$X
y = data$Y
theta.hat = para.solve(x, y, Model = 1)
y.tilde = y.tilde.sample(x, theta.hat, Model = 1)
Tn.hat = Tnsc(x, y, y.tilde, c = 5)
}
