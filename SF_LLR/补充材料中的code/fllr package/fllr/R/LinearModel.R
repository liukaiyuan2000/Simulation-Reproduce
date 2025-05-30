#' Functional linear model
#'
#' Estimates the standard scalar-on-function linear model, where all functional
#' data are projected into a subspace given by basis functions.
#'
#' @param Responses vector of scalar responses
#' @param LEARN     matrix of regressors, functional values evaluated at common
#' domain, one function per row
#' @param PRED      matrix of functions at which new values are predicted, one
#' function per row
#' @param percent   if \code{BASIS} is not given, percentage of explained sample
#' variance that is taken into account when building the basis expansion in
#' terms of eigenfunctions of the empirical covariance operator
#' @param BASIS     if basis functions are not estimated using the empirical
#' covariance operator, matrix of basis functions. One function per row.
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{Predicted.responses} A vector of predicted responses, each
#' element of the vector corresponding to a row of \code{PRED}.
#' \item \code{PREDICTED.DERIV} A matrix of predicted (Riesz representations of)
#' functional derivatives, each row corresponding to a row of \code{PRED}.
#' \item \code{mse} Minimum mean squared error obtained in the search
#' leave-one-out cross-validation for the optimal dimension of the basis for
#' the regression operator.
#' \item \code{dim.opt} Optimal basis dimension used for both the regression
#' operator and the functional derivative. If \code{BASIS} is provided,
#' \code{dim.opt} is chosen as the numer of rows of \code{BASIS}.
#' }
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#' @references Reiss, P. T., and Ogden, R. T. (2007). Functional principal
#' component regression and functional partial least squares. \emph{Journal of
#' the American Statistical Association}. 102, 984--996.
#'
#' @examples
#' dat = generate(100,200,4)
#' LEARN = dat$LEARN
#' PRED = dat$PRED
#' Responses = dat$Responses
#'
#' res = LinearModel(Responses, LEARN, PRED)
#'
#' plot(res$Predicted.responses~dat$oracle.Responses)
#'

LinearModel = function(Responses, LEARN, PRED, percent = 0.95, BASIS = NULL)
{
  outsample = T
  if(identical(LEARN,PRED)) outsample = F
  ## if LEARN is the same matrix as PRED, outsample = F and PRED = LEARN
  ##
  # Estimate the second order moment operator
  nlearn = nrow(LEARN)
  gridsize = ncol(LEARN)
  if(is.null(BASIS)){
    if(outsample){
      sample.size = nrow(LEARN) + nrow(PRED)
      # COV = crossprod(rbind(LEARN, PRED)) / (sample.size * (gridsize - 1))
      COV = var(rbind(LEARN,PRED)) # added in R1
    }else{
      sample.size = nrow(LEARN)
      # COV = crossprod(LEARN) / (sample.size * (gridsize - 1))
      COV = var(LEARN) # added in R1
    }
    # Compute estimated first "dim.max" eigenfunctions
    res.eig = eigen(COV, sym = T)
    # Drop first eigenvalue and derive the approximating subspace expressing at
    # least "percent" of the remaining total variance
    cum.part.of.variance = cumsum(res.eig$values[-1]) / sum(res.eig$values[-1])
    dim.max = sum(cum.part.of.variance < percent) + 2
    BASIS = t(res.eig$vectors[, 1:dim.max]) * sqrt(gridsize - 1)
  } else { dim.max = nrow(BASIS) }

  INNER.LEARN = tcrossprod(LEARN, BASIS) / (gridsize - 1)
  H = matrix(Inf, nrow=nrow(INNER.LEARN),ncol=1)
  MSE = llsm_cv(matrix(INNER.LEARN,ncol=dim.max),matrix(0,nrow=nlearn,ncol=
                                 nlearn),Responses,H=H,kernel="Epanechnikov")$CV
  dim.opt = order(MSE)[1]-1
  mse = min(MSE)
  nout = nrow(PRED)

  if(dim.opt==0){
    Predicted.responses = rep(mean(Responses),nout)
  # llsm(NULL,NULL,matrix(0,nrow=nout,ncol=nlearn),Responses,h=Inf,kernel="Epa")
    PREDICTED.DERIV = matrix(0,nrow=nout,ncol=gridsize)
  } else {
    B = matrix(BASIS[1:dim.opt,],nrow=dim.opt)
    INNER.LEARN = matrix(INNER.LEARN[,1:dim.opt],ncol=dim.opt)
    INNER.PRED = tcrossprod(PRED, B) / (gridsize - 1)
    # res = llsm(INNER.LEARN,INNER.PRED,matrix(0,nrow=nout,ncol=nlearn),
    # Responses,h=Inf,kernel="Epa")
    # summary(cbind(1,INNER.PRED)%*%linmod$coef-res[,1])
    linmod = lm(Responses~INNER.LEARN)
    Predicted.responses = cbind(1,INNER.PRED)%*%linmod$coef
    PREDICTED.DERIV = matrix(rep(linmod$coef[-1],nout),nrow=nout,byrow=TRUE)%*%B
  }
  return(list(Predicted.responses = Predicted.responses,
              PREDICTED.DERIV = PREDICTED.DERIV,
              mse = mse,
              dim.opt = dim.opt))
}
