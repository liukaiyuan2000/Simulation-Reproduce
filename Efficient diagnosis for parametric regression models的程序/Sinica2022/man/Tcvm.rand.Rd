% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Funcs.R
\name{Tcvm.rand}
\alias{Tcvm.rand}
\title{Calculate Tcvm.rand}
\usage{
Tcvm.rand(V, theta, e)
}
\arguments{
\item{V}{cbind(X.hat, Z) a n*(p+q) matrix}

\item{e}{e.hat}
}
\value{
Tcvm test statistics(2.5)
}
\description{
Calculate Tcvm.rand
}
\examples{
Tcvm.rand(V, theta.normal, e.hat)
}
