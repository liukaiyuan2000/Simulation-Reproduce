#' Gaussian process regression
#'
#' Estimates the scalar-on-function linear model using a Gaussian process method
#' for a single functional predictor, and a scalar response.
#'
#' @param Responses vector of scalar responses
#' @param LEARN     matrix of regressors, functional values evaluated at common
#' domain, one function per row
#' @param PRED      matrix of functions at which new values are predicted, one
#' function per row
#'
#' @return \code{Predicted.responses} A vector of predicted responses, each
#' element of the vector corresponding to a row of \code{PRED}.
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#' @references Wang, B., and Xu, A. (2019).
#'	Gaussian process methods for nonparametric functional regression with mixed
#'	predictors.	\emph{Computational Statistics & Data Analysis}, 131, 80--90.
#'
#' @examples
#' dat = generate(100,200,4)
#' LEARN = dat$LEARN
#' PRED = dat$PRED
#' Responses = dat$Responses
#'
#' res = GPR(Responses, LEARN, PRED)
#'
#' plot(res$Predicted.responses~dat$oracle.Responses)
#'

GPR = function(Responses, LEARN, PRED){
  # Gaussian process regression with functional predictor
  # method taken from Wang and Xu (CSDA, 2019)
  gridsize = ncol(LEARN)
  require(parallelDist)
  DIST_full = as.matrix(dist(rbind(LEARN,PRED), method = "euclidean"))/
    sqrt(gridsize-1)
  DIST = DIST_full[1:nrow(LEARN),1:nrow(LEARN)]
  # Covariance function of teh Gaussian process
  K = function(v,eta,sigma){
    v^2*exp(-eta^2*DIST) + sigma^2*diag(1,nrow=nrow(DIST))
  }
  num_inv = function(A){
    # numerical inverse of a matrix that does not produce errors
    Ai = try(solve(A),silent=TRUE)
    if(class(Ai)=="try-error"){
      warning("Numerically singular matrix in GPR, using inifite matrix
              instead")
      return(diag(Inf,nrow(A)))
    } else return(Ai)
  }
  # log-likelihood function to maximize
  lik = function(par){
    c(determinant(K(par[1],par[2],par[3]),log=TRUE)$modulus+
        Responses%*%num_inv(K(par[1],par[2],par[3]))%*%Responses)
  }
  pars = suppressWarnings(optim(c(1,1,1),lik)$par)
  Ki = num_inv(K(pars[1],pars[2],pars[3]))
  DIST2 = DIST_full[-(1:nrow(LEARN)),1:nrow(LEARN)]
  # in rows are distances of elements of PRED from all LEARN
  v = pars[1];  eta = pars[2];  sigma = pars[3];
  Ks = v^2*exp(-eta^2*DIST2)
  pe = c(Ks%*%Ki%*%Responses)
  return(list(Predicted.responses=pe))
}
