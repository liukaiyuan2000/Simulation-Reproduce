#' Eigenbasis selection
#'
#' Automated eigenbasis selection for densely observed functional data.
#'
#' @param X data matrix, functional values evaluated at common domain,
#' one function per row
#' @param method  method used for the selection of the number of eigenvectors,
#' possible values are \code{FVE} for the fraction-of-variance-explained in the
#' fitted covariance function, or \code{BIC} for the conditional pseudo-Gaussian BIC selection
#' @param FVEthreshold  a numerical constant between 0 and 1, giving the threshold
#' for method \code{FVE}. Ignored with \code{method=BIC}.
#'
#' @return A matrix whose columns correspond to the estimated basis of functions.
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#' @examples
#' gridsize = 101                           # grid size
#' Grid = seq(0, 1, length = gridsize)      # grid of measurements
#' n = 150                                  # sample size
#' Jmodel = 4                               # number of basis elements used for building functional predictors
#' BASIS = phif(Jmodel, Grid)               # first Jmodel Fourier basis elements
#' UNIF = matrix(runif(n * Jmodel, min = -1, max = 1), nrow = n, ncol = Jmodel)
#' X = UNIF %*% BASIS[1:Jmodel, ]
#'
#' # dimensions selected using different methods
#' dim(eigenbasis(X))
#' dim(eigenbasis(X,FVEthreshold = .95))
#' dim(eigenbasis(X,method="BIC"))

eigenbasis = function(X, method="FVE", FVEthreshold=.9){
  # finds an optimal FPCA basis for densely observed functional data
  # based on
  #    FVE: fraction-of-variance-explained threshold the fitted covariance, or
  #    BIC: conditional pseudo-Gaussian BIC selection of J for FPCA
  method = match.arg(method,c("BIC","FVE"))
  mu = colMeans(X)          # estimated mean function for dense functional data
  Xc = t(t(X)-mu)           # centred functional data
  V = var(X)
  # all.equal(var(Y),crossprod(Yc)/(nrow(Y)-1))
  eX = eigen(V)
  B = t(eX$vectors)
  if(method=="BIC"){
    sigma2 = mean(diag(V))    # estimate of the residual variance
    IC = rep(NA,ncol(X))
    C = tcrossprod(Xc,B)
    N = nrow(X)*ncol(X)
    for(K in 1:ncol(X)){
      Xhat = t(t(C[,1:K,drop=FALSE]%*%B[1:K,,drop=FALSE])+mu)
      IC[K] = sum(diag(tcrossprod(X-Xhat)))/(sigma2) + N*log(sigma2)+K*log(N)
      if(K > 1 && IC[K] > IC[K-1]){
        # cease whenever BIC stops decreasing as in PACE
        oK = K-1;
        break;
      } else if(K == ncol(X)){
        oK = K
      }
    }
  }
  if(method=="FVE"){
    cumFVE = cumsum(eX$values)/sum(eX$values)
    oK = min(which(cumFVE > FVEthreshold))
  }
  return(t(B[1:oK,,drop=FALSE]))
}
