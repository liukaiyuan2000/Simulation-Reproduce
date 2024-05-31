#' @title Data Generate Process
#'
#' @param N sample size N
#' @param m dimensions
#' @param case simulation cases, Table 1-4
#'
#' @return   a N*m matrix
#' @export
#'
#' @examples DGP(N = 4, m = 4, case = 1)
DGP <- function(N, m, case){
  Sigma <- matrix(0.1, m, m)
  diag(Sigma) <- 1
  X <- switch(
    case, 
    matrix(rnorm(N * m), N, m), 
    matrix(rcauchy(N * m), N, m), 
    matrix(rt(N * m, 2), N, m), 
    mvrnorm_rcpp(N, rep(0, m), Sigma)
  )
  X
}