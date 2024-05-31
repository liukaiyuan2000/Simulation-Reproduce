#' @title Data Generate Process
#'
#' @param n sample size N = n + 1
#' @param m dimensions
#' @param rho correlation coefficients when i not equal to j
#'
#' @return   a n*m matrix
#' @export
#'
#' @examples DGP(n = 5, m = 4, rho = 0.1)
DGP <- function(n, m, rho){
  if(rho == 0){
    X <- matrix(rnorm(n*m), n, m)
  }else{
    Sigma <- matrix(rho, m, m)
    diag(Sigma) <- 1
    X <- mvrnorm_rcpp(n, rep(0, m), Sigma)
  }
  return(X)
}
