#' Local constant functional estimation
#'
#' Estimates the scalar-on-function model using the local constant
#' (Nadaraya-Watson's estimator) for the regression operator and functional
#' data. Does not directly provide an estimator of functional derivatives.
#' Bandwidth is chosen via leave-one-out cross-validation.
#'
#' @param Responses vector of scalar responses
#' @param DIST      a square symmetric matrix of all distances of regressor
#' functions from each other, see \code{\link{L2metric}}.
#' @param DIST.PRED an \code{m}-times-\code{n} matrix of distances of all
#' functions at which we intend to predict (\code{m} functions) from the
#' regressor functions (\code{n} functions).
#' @param kNN.grid.length  grid length for the bandwidth search.
#' @param fullgrid  logical; should the full grid [0,1] of k-NN be used
#' directly for cross-validation? If \code{FALSE}, the grid is built iteratively
#' from small to larger \code{k}. By default, \code{fullgrid=FALSE}.
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{Predicted.responses} A vector of predicted responses, each
#' element of the vector corresponding to a function that corresponds to a row
#' in \code{DIST}.
#' \item \code{mse} Minimum mean squared error obtained in the search for the
#' optimal bandwidth for the regression operator.
#' \item \code{knn.opt} The optimal bandwidths for the regression operator
#' obtained by leave-one-out cross-validation.
#' }
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#' @references Ferraty, F., and Vieu, P. (2006). Nonparametric functional data
#' analysis. \emph{Springer Series in Statistics. Springer, New York}.
#'
#' @examples
#' dat = generate(100,200,4)
#' LEARN = dat$LEARN
#' PRED = dat$PRED
#' Responses = dat$Responses
#'
#' DIST = L2metric(LEARN,LEARN)
#' DIST.PRED = L2metric(PRED,LEARN)
#' res = KernelSmoother(Responses, DIST, DIST.PRED)
#'
#' plot(res$Predicted.responses~dat$oracle.Responses)
#'

KernelSmoother = function(Responses, DIST, DIST.PRED, kNN.grid.length = 30,
                          fullgrid=FALSE)
{
  nlearn = nrow(DIST)
  DIST.SORTED.ROW.BY.ROW = t(apply(DIST, 1, sort))
  DIST.PRED.SORTED.ROW.BY.ROW = t(apply(DIST.PRED, 1, sort))
  #
  if(nlearn<50) warning("With less than 50 learning curves the cross-validation
                        for the kNN method may be unstable.")
  if(fullgrid) Knn.learn.ratio = c(0.02,1) else
    Knn.learn.ratio = c(0.02, 0.25, 0.5, 0.75, 1)
  # Build an adaptative grid of kNN (number of knearest neighbours used for
  # determining local bandwidths). If the minimum of the criterion ("CRIT.REG")
  # is reached for the bandwidth corresponding to "Knn.grid[kNN.grid.length]",
  # then the grid "Knn.grid" is updated
  kNN.iterate = TRUE # indicator stopping the iteration in kNN CV procedure
  knn.max = 2        # initiate the kNN CV procedure
  ratio = 1
  while(kNN.iterate){
    knn.min = min(max(round(nlearn * Knn.learn.ratio[ratio]), knn.max),nlearn-1)
    knn.max = min(max(round(nlearn * Knn.learn.ratio[ratio + 1]),
                      knn.min + kNN.grid.length  - 1),nlearn-1)
    ratio = ratio + 1
    if(ratio == length(Knn.learn.ratio)-1) kNN.iterate = FALSE # stop iteration
    if(knn.max == nlearn-1) kNN.iterate = FALSE                # stop iteration
    knn.range = knn.max - knn.min
    if(knn.range <= kNN.grid.length){
      # this now happens only if knn.max == nlearn
      Knn.grid = knn.min:knn.max
    }else{
      step = ceiling(knn.range / kNN.grid.length)
      Knn.grid = seq(from = knn.min, to = knn.max, by = step)
    }
    actual.kNN.grid.length = length(Knn.grid)
    rank.min = actual.kNN.grid.length
    # Knn.learn.ratio = c(0.02, 0.25, 0.5, 0.75, 0.9)
    # Build an adaptive grid of kNN (number of k nearest neighbors used for
    # determining local bandwidths). If the minimum of the criterion ("MSE")
    # is reached for the bandwidth corresponding to "Knn.grid[kNN.grid.length]",
    # then the grid "Knn.grid" is updated
    ############################################################################
    # estimate regression operator with bandwidth maximizing the CV criterion
    ############################################################################
    MSE = rep(0, kNN.grid.length)
    H = matrix(nrow=nlearn,ncol=length(Knn.grid))
    for(knni in 1:length(Knn.grid)) H[,knni] = DIST.SORTED.ROW.BY.ROW[,
                                                            Knn.grid[knni] + 1]
    MSE = ksm_cv(DIST,Responses,H=H, kernel="Epanechnikov")$CV
    if(rank.min != order(MSE)[1]) break
    # if minimum CV is not the last one, break
  }
  # Determine corresponding criterion value
  mse = min(MSE)
  # Determine optimal number of k-nearest neighbours and dimension J
  knn = (Knn.grid)[order(MSE)[1]]
  #
  Predicted.responses = ksm(DIST.PRED, Responses, DIST.PRED.SORTED.ROW.BY.ROW[,
                                                                      knn + 1])
  return(list(Predicted.responses = Predicted.responses,
              mse = mse,
              knn.opt = knn))
}

#' Cross-validation for the local constant smoother
#'
#' For the local constant smoother (Nadaraya-Watson's estimator)
#' in the scalar-on-function linear model, given a sequence of bandwidths,
#' returns a sequence of (leave-one-out) cross-validated mean squared errors
#' for bandwidth selection.
#'
#' @param D         symmetric square matrix of distances of the functional data
#' from each other, one row per function
#' @param Y         vector of scalar responses
#' @param H         bandwidths to be considered. Could be either a vector of
#' bandwidths that will be applied to each function, or a matrix with individual
#' bandwidths applied to learning function, one row per function. Alternatively,
#' if \code{kNN} is specified, the bandwidth matrix is built automatically using
#' the local nearest neighbors bandwidths.
#' @param kNN       sequence of nearest neighbors (as a fraction of total number
#' of neighbors, i.e. a number between 0 and 1) to be considered for local
#' bandwidth selection. If \code{kNN} is provided, matrix \code{H} is build
#' automatically by evaluating sample quantiles of the distances in rows of
#' the matrix of distances \code{D}.
#' @param nCV       if provided, only the first \code{nCV} leave-one-out
#' cross-validated trials is considered in the resulting mean squared error.
#' By default, all functions are included in cross-validation.
#' @param kernel    kernel function. Possible choices are \code{"uniform"},
#' \code{"triangular"}, \code{"Epanechnikov"}, \code{"biweight"},
#' \code{"triweight"}, and \code{"Gaussian"}.
#' @param quantile.type if \code{kNN} is provided, parameter that enters
#' function \code{\link{quantile}} that specifies the type of the sample
#' quantile taken.
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{CV} A matrix that corresponds to \code{H} that contains all the
#' resulting mean squared errors of the cross-validation procedure.
#' \item \code{kernel} Kernel used.
#' }
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#' @references Ferraty, F., and Vieu, P. (2006). Nonparametric functional data
#' analysis. \emph{Springer Series in Statistics. Springer, New York}.
#'
#' @examples
#' dat = generate(100,200,4)
#' LEARN = dat$LEARN
#' PRED = dat$PRED
#' Responses = dat$Responses
#'
#' ksm_cv(L2metric(LEARN, LEARN), Responses, kNN = seq(.1,.9,length=11),
#' kernel = "Gauss")
#'

ksm_cv <- function(D,Y,H=NULL,kNN=NULL,nCV=NULL, kernel=c("uniform",
  "triangular","Epanechnikov","biweight","triweight","Gaussian"),
  quantile.type=7) {
  krn = c("uniform","triangular","Epanechnikov","biweight","triweight",
          "Gaussian")
  kernel = match.arg(kernel)
  kernI = which(kernel==krn)
  if(is.null(H) & is.null(kNN)) stop("Neither candidate bandwidths nor nearest
                                     neighbour proportions specified.")
  if((!is.null(H)) & (!is.null(kNN))) stop("Specify either bandwidths or nearest
                                           neighbour proportions, not both.")
  # construct full matrix of candidate bandwidths, different for each function
  if(is.null(H)) Heval = length(kNN)
  if(is.null(kNN)){ if(is.null(dim(H))) Heval = length(H) else Heval = ncol(H) }
  Hf = matrix(nrow=nrow(D),ncol=Heval)
  # each function has now its own sequence of candidate bw
  for(i in 1:nrow(D)){
    if(is.null(H)){
      Hf[i,] = quantile(D[i,-i],kNN,type=quantile.type)
      Hf[i,kNN==1] = Inf
    }
    if(is.null(kNN)){
      if(is.null(dim(H))) Hf[i,] = H else Hf[i,] = H[i,]
    }
  }
  Hf[!is.finite(Hf)] = -1
  if(is.null(nCV)) nCV = nrow(D)
  if(nCV>nrow(D)) nCV = nrow(D)
  out <- .C("ksm_cv", D = as.double(D), Y = as.double(Y),
            n = as.integer(nrow(D)),
            Hf = as.double(Hf), nH = as.integer(Heval),
            CV=as.double(rep(0,Heval)), nCV = as.integer(nCV),
            kernI = as.integer(kernI))
  CV = out$CV
  Hf[Hf==-1] = Inf # back to the usual coding for Inf
  if(!is.null(H)){ # H cross-validation
    return(list(CV = CV, kernel = kernel))
  }
  if(is.null(H)){ # kNN cross-validation
    return(list(CV = CV, kernel = kernel))
  }
}

#' Local constant smoother with a fixed bandwidth
#'
#' The local constant smoother (Nadaraya-Watson's estimator)
#' in the scalar-on-function linear model with a given bandwidth.
#'
#' @param D         a (possibly non-symmetric) matrix of distances of the
#' learning sample of functional data from the predictor functions, one row per
#' a learning sample function
#' @param Y         vector of scalar responses
#' @param h         a single bandwidth. Alternatively, if \code{kNN} is
#' specified, local bandwidths are built automatically using the nearest
#' neighbors.
#' @param kNN       number of nearest neighbors (as a fraction of total number
#' of neighbors, i.e. a number between 0 and 1) to be considered for local
#' bandwidth selection. If \code{kNN} is provided, vector \code{h} is build
#' automatically by evaluating sample quantiles of the distances in rows of the
#' matrix of distances \code{D}.
#' @param kernel    kernel function. Possible choices are \code{"uniform"},
#' \code{"triangular"}, \code{"Epanechnikov"}, \code{"biweight"},
#' \code{"triweight"}, and \code{"Gaussian"}.
#' @param quantile.type if \code{kNN} is provided, parameter that enters
#' function \code{\link{quantile}} that specifies the type of the sample
#' quantile taken.
#'
#' @return A vector of predicted responses whose length is the number of columns
#' of \code{D}.
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#' @references Ferraty, F., and Vieu, P. (2006). Nonparametric functional data
#' analysis. \emph{Springer Series in Statistics. Springer, New York}.
#'
#' @examples
#' dat = generate(100,200,4)
#' LEARN = dat$LEARN
#' PRED = dat$PRED
#' Responses = dat$Responses
#'
#' ksm(L2metric(LEARN, PRED), Responses, kNN = .2, kernel = "Gauss")
#'

ksm <- function(D, Y, h=NULL, kNN = NULL, kernel=c("uniform","triangular",
           "Epanechnikov","biweight","triweight","Gaussian"), quantile.type=7) {

  if(is.null(dim(D))) D = matrix(D,nrow=1)
  if(is.null(h) & is.null(kNN)) stop("Neither bandwidth nor nearest neighbour
                                     proportion specified.")
  if((!is.null(h)) & (!is.null(kNN))) stop("Specify either bandwidth or nearest
                                           neighbour proportion, not both.")
  if(is.null(h) & (!is.null(kNN))){
    h = apply(D,1,function(x) quantile(x,kNN,type=quantile.type))
    if(kNN == 1) h = rep(Inf,nrow(D))
  }
  if(length(h)==1) h = rep(h,nrow(D))
  h[!is.finite(h)] = -1
  krn = c("uniform","triangular","Epanechnikov","biweight","triweight",
          "Gaussian")
  kernel = match.arg(kernel)
  kernI = which(kernel==krn)
  out <- .C("ksm", D = as.double(D), Y = as.double(Y), h = as.double(h),
            n = as.integer(length(Y)), m = as.integer(nrow(D)),
            kernI = as.integer(kernI), res=as.double(rep(0,nrow(D))))
  return(out$res)
}
