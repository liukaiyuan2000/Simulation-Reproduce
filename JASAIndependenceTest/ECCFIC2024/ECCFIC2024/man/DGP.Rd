% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DGP.R
\name{DGP}
\alias{DGP}
\title{Data generate process}
\usage{
DGP(
  example,
  model = NULL,
  n = NULL,
  p = NULL,
  q = NULL,
  delta = NULL,
  slice.num = NULL
)
}
\arguments{
\item{example}{example 1-10}

\item{model}{choose model, for example 1 and 7}

\item{n}{sample size}

\item{p}{dimension of \eqn{X}}

\item{q}{dimension of \eqn{Y}}

\item{delta}{parameter \eqn{\delta} in example 2}

\item{slice.num}{number of slices in example 4}
}
\value{
a list contain \eqn{X, Y}, and the discrete \eqn{Y}.
}
\description{
Data generate process
}
\examples{
example = 1
model = 1
n = 25
p = 5
q = 1
dat = DGP(exmaple = example, model = model, n = 25, p = 5, q = 1)
}
