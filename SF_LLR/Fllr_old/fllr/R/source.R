#' Oracle bandwidth value for functional derivative estimation
#'
#' Given the true (Riesz representatives of) functional derivatives in a scalar-on-function
#' regression problem, this function evaluates the oracle bandwidth choice for the estimation
#' of the functional derivatives in a given functional basis.
#'
#' @param Responses vector of scalar responses
#' @param LEARN     matrix of regressors, functional values evaluated at common domain,
#' one function per row
#' @param ORACLE.DERIV matrix of the true functional derivatives, each row corresponding
#' to a row of \code{LEARN}.
#' @param BASIS     basis functions into which the random functions are decomposed. One
#' function per row. The true dimension of this space is assumed to be the number of rows
#' of this matrix.
#' @param kNN.grid.length  grid length for the bandwidth search
#'
#' @return A single bandwidth in terms of the number of nearest neighbors that
#' corresponds to \code{knn.opt$deriv} from \code{\link{fllr}}. Works for
#' \code{method = cv} only.
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#'	@references Ferraty, F., and Nagy, S. (2021).
#'	Scalar-on-function local linear regression and beyond.
#'	\emph{Biometrika}, to appear.

oracle.bw = function(Responses, LEARN, ORACLE.DERIV, BASIS, kNN.grid.length = 30,
                     fullgrid=FALSE)
{
  if(is.null(fullgrid)) fullgrid = (nrow(LEARN)<=250)
  gridsize = ncol(LEARN)
  DIST = as.matrix(dist(LEARN, method = "euclidean")) / sqrt(gridsize - 1)
  DIST.SORTED.ROW.BY.ROW = t(apply(DIST, 1, sort))
  nlearn = nrow(LEARN)
  B = BASIS
  dim = nrow(BASIS)
  INNER.LEARN = tcrossprod(LEARN, B) / (gridsize - 1)
  PILOT.BOOT.DERIV = ORACLE.DERIV
  # Determine the optimal bandwidth for estimating the functional derivatives
  if(nlearn<50) warning("With less than 50 learning curves the cross-validation for the kNN method may be unstable.")
  if(fullgrid) Knn.learn.ratio = c(0.02,1) else
    Knn.learn.ratio = c(0.02, 0.25, 0.5, 0.75, 1)
  # Build an adaptative grid of kNN (number of knearest neighbours used for determining local bandwidths)
  # If the minimum of the criterion ("CRIT.REG") is reached for the bandwidth corresponding to
  # "Knn.grid[kNN.grid.length]", then the grid "Knn.grid" is updated
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
    #########################################################################################
    # estimate kNN
    #########################################################################################
    Crit.deriv.boot = rep(NA,kNN.grid.length)
    count = 0
    for(knn in Knn.grid){
      count = count + 1
      COEF.DERIV = llsm_leave(INNER.LEARN,INNER.LEARN,DIST,Responses,h=DIST.SORTED.ROW.BY.ROW[, knn + 1],kernel="Epa")[,-1]
      ESTIMATED.DERIV = COEF.DERIV %*% B
      # Compute the criterion for estimating the fuctional derivative
      Crit.deriv.boot[count] = mean((PILOT.BOOT.DERIV - ESTIMATED.DERIV)^2)
    }
    if(rank.min != order(Crit.deriv.boot)[1]) break
  }
  # Determine the criterion value
  crit.deriv.boot = min(Crit.deriv.boot)
  # Determine the optimal number of k-nearest neighbours for the functional derivative
  knn.deriv = (Knn.grid)[which.min(Crit.deriv.boot)]
  # Save optimal bandwidths in terms of k-nearest neighbours into outputs
  knn.opt = knn.deriv
  return(knn.opt)
}

###
# fllr
###

matsolve <- function(a, b)
{
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(b))
  stopifnot(is.finite(a))
  stopifnot(is.finite(b))
  stopifnot(is.matrix(a))
  # stopifnot(a == t(a))
  b <- as.matrix(b)
  stopifnot(nrow(a) == nrow(b))
  storage.mode(a) <- "double"
  storage.mode(b) <- "double"
  .C("matsolve", a = a, b = b, nrowb = nrow(b), ncolb = ncol(b),
     PACKAGE = "fllr")$b
}

matvecmult <- function(a, x)
{
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(x))
  stopifnot(is.finite(a))
  stopifnot(is.finite(x))
  stopifnot(is.matrix(a))
  stopifnot(ncol(a) == length(x))
  .C("matvecmult", a = as.double(a), x = as.double(x),
     nrow = nrow(a), ncol = ncol(a), result = double(nrow(a)),
     PACKAGE = "fllr")$result
}

matmatmult <- function(a, b)
{
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(b))
  stopifnot(is.finite(a))
  stopifnot(is.finite(b))
  stopifnot(is.matrix(a))
  stopifnot(is.matrix(b))
  stopifnot(ncol(a) == nrow(b))
  .C("matmatmult", a = as.double(a), b = as.double(b),
     nrowa = nrow(a), ncola = ncol(a), ncolb = ncol(b),
     c = matrix(as.double(0), nrow = nrow(a), ncol = ncol(b)),
     PACKAGE = "fllr")$c
}

#' Functional local linear smoother with a fixed bandwidth
#'
#' The functional local linear smoother
#' in the scalar-on-function linear model with a given bandwidth.
#'
#' @param C         a matrix of the basis coefficients of the the learning sample
#' functions (n functions), one row per function
#' @param Cnew      a matrix of the basis coefficients of the predictor functions
#' (m functions), one row per predictor
#' @param D         a matrix of distances of the predictor functions from the
#' learning sample functions (m-times-n matrix), one row per a predictor
#' @param Y         vector of scalar responses
#' @param h         either a single bandwidth, or a vector of m bandwidths to be used
#' for the predictor functions. Alternatively, if \code{kNN} is specified,
#' local bandwidths are built automatically using the nearest neighbors.
#' @param kNN       number of nearest neighbors (as a fraction of total number of neighbors, i.e.
#' a number between 0 and 1) to be considered for local bandwidth
#' selection. If \code{kNN} is provided, vector \code{h} is build automatically by
#' evaluating sample quantiles of the distances in rows of the matrix of distances \code{D}.
#' @param kernel    kernel function. Possible choices are \code{"uniform"}, \code{"triangular"},
#' \code{"Epanechnikov"}, \code{"biweight"}, \code{"triweight"}, and \code{"Gaussian"}.
#' @param quantile.type if \code{kNN} is provided, parameter that enters function \code{\link{quantile}}
#' that specifies the type of the sample quantile taken
#'
#' @return A matrix where each row corresponds to a single predictor.
#' The columns are the estimated regression coefficients in the given basis expansion.
#' The first column is the intercept - an estimate of the regression operator. Remaining
#' columns correspond to the basis expansion of the estimator of the functional derivative
#' of the regression operator.
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#'	@references Ferraty, F., and Nagy, S. (2021).
#'	Scalar-on-function local linear regression and beyond.
#'	\emph{Biometrika}, to appear.
#'
#' @examples
#' dat = generate(100,200,4)
#' LEARN = dat$LEARN
#' PRED = dat$PRED
#' Responses = dat$Responses
#' BASIS = dat$BASIS
#'
#' C = project(LEARN,BASIS)
#' Cnew = project(PRED,BASIS)
#' D = L2metric(PRED,LEARN)
#' res = llsm(C, Cnew, D, Responses, kNN = .05, kernel = "Gauss")
#'
#' plot(res[,1]~Regression[-(1:nlearn)])

llsm <- function(C, Cnew, D, Y, h=NULL, kNN = NULL, kernel=c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian"), quantile.type=7) {
  if(is.null(dim(D))){
    D = matrix(D,nrow=1)
  }
  krn = c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian")
  kernel = match.arg(kernel,krn)
  kernI = which(kernel==krn)
  if(is.null(h) & is.null(kNN)) stop("Neither bandwidth nor nearest neighbour proportion specified.")
  if((!is.null(h)) & (!is.null(kNN))) stop("Specify either bandwidth or nearest neighbour proportion, not both.")
  if(is.null(h) & (!is.null(kNN))){
    h = apply(D,1,function(x) quantile(x,kNN,type=quantile.type))
    if(kNN == 1) h = rep(Inf,nrow(D)) # if all functions should be used perform usual (linear) regression
  }
  # if(h==Inf) return(llsm(C,Cnew,rep(1,length(D)),Y,h=-1,allJ = allJ, kernel = kernel)) # if h==Inf return usual regression output
  h[!is.finite(h)] = -1
  if(is.null(C)){
    iC = matrix(1,nrow=ncol(D))
    iCnew = matrix(1,nrow=nrow(D),1)
  } else {
    iC<-cbind(1,C)
    iCnew = matrix(NA,nrow(D),ncol(iC))
    iCnew[,1] = 1
    iCnew[,-1] = Cnew
  }
  out <- .C("llsm_single", C = as.double(iC), Cnew = as.double(iCnew), D = as.double(D), Y = as.double(Y),
            h = as.double(h),
            n = as.integer(nrow(iC)), m = as.integer(nrow(D)), J = as.integer(ncol(iC)),
            kernI = as.integer(kernI),
            res=as.double(rep(0,nrow(D)*ncol(iC))))
  return(matrix(out$res,nrow=nrow(D)))
}

lqsm <- function(C, Cnew, D, Y, h=NULL, kNN = NULL, kernel=c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian"), quantile.type=7) {
  if(is.null(dim(D))){
    D = matrix(D,nrow=1)
  }
  krn = c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian")
  kernel = match.arg(kernel,krn)
  kernI = which(kernel==krn)
  if(is.null(h) & is.null(kNN)) stop("Neither bandwidth nor nearest neighbour proportion specified.")
  if((!is.null(h)) & (!is.null(kNN))) stop("Specify either bandwidth or nearest neighbour proportion, not both.")
  if(is.null(h) & (!is.null(kNN))){
    h = apply(D,1,function(x) quantile(x,kNN,type=quantile.type))
    if(kNN == 1) h = rep(Inf,nrow(D)) # if all functions should be used perform usual (linear) regression
  }
  # if(h==Inf) return(llsm(C,Cnew,rep(1,length(D)),Y,h=-1,allJ = allJ, kernel = kernel)) # if h==Inf return usual regression output
  h[!is.finite(h)] = -1
  if(is.null(C)){
    iC = matrix(1,nrow=ncol(D))
    iCnew = matrix(1,nrow=nrow(D),1)
  } else {
    iC<-cbind(1,C)
    iCnew = matrix(NA,nrow(D),ncol(iC))
    iCnew[,1] = 1
    iCnew[,-1] = Cnew
  }
  out <- .C("lqsm_single", C = as.double(iC), Cnew = as.double(iCnew), D = as.double(D),
            Y = as.double(Y),
            h = as.double(h),
            n = as.integer(nrow(iC)), m = as.integer(nrow(D)),
            kernI = as.integer(kernI),
            res=as.double(rep(0,nrow(D)*3)))
  return(matrix(out$res,nrow=nrow(D)))
}

llsm_leave <- function(C, Cnew, D, Y, h=NULL, kNN = NULL, kernel=c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian"), quantile.type=7) {
  ## input:
  # C:    basis coefficients of the X functions (matrix n*J)
  # Cnew: basis coefficients of the new regressor x0 (matrix m*J)
  # D:    matrix of distances of x0 from the X functions (matrix m*n)
  # Y:    response vector (vector n)
  # h:    bandwidth (vector m or vector 1, in latter case replicated)
  # kNN:  proportion of k-nearest neighbors to be used in bandwidth definition (in the interval (0,1))

  if(is.null(dim(D))){
    D = matrix(D,nrow=1)
  }
  krn = c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian")
  kernel = match.arg(kernel)
  kernI = which(kernel==krn)
  if(is.null(h) & is.null(kNN)) stop("Neither bandwidth nor nearest neighbour proportion specified.")
  if((!is.null(h)) & (!is.null(kNN))) stop("Specify either bandwidth or nearest neighbour proportion, not both.")
  if(is.null(h) & (!is.null(kNN))){
    h = apply(D,1,function(x) quantile(x,kNN,type=quantile.type))
    if(kNN == 1) h = rep(Inf,nrow(D)) # if all functions should be used perform usual (linear) regression
  }
  # if(h==Inf) return(llsm(C,Cnew,rep(1,length(D)),Y,h=-1,allJ = allJ, kernel = kernel)) # if h==Inf return usual regression output
  h[!is.finite(h)] = -1
  if(is.null(C)){
    iC = matrix(1,nrow=ncol(D))
    iCnew = matrix(1,nrow=nrow(D),1)
  } else {
    iC<-cbind(1,C)
    iCnew = matrix(NA,nrow(D),ncol(iC))
    iCnew[,1] = 1
    iCnew[,-1] = Cnew
  }
  out <- .C("llsm_single_leave", C = as.double(iC), Cnew = as.double(iCnew), D = as.double(D), Y = as.double(Y),
            h = as.double(h),
            n = as.integer(nrow(iC)), m = as.integer(nrow(D)), J = as.integer(ncol(iC)),
            kernI = as.integer(kernI),
            res=as.double(rep(0,nrow(D)*ncol(iC))))
  return(matrix(out$res,nrow=nrow(D)))
}

#' Cross-validation for the functional local linear smoother
#'
#' For the functional local linear smoother
#' in the scalar-on-function linear model, given a sequence of bandwidths,
#' returns a sequence of (leave-one-out) cross-validated mean squared errors
#' for bandwidth selection.
#'
#' @param C         a matrix of the basis coefficients of the the learning sample
#' functions (n functions), one row per function
#' @param D         a symmetric square matrix of distances of the learning sample
#' functions from themselves (n-times-n matrix), one row per a function
#' @param Y         vector of scalar responses
#' @param H         bandwidths to be considered. Could be either a vector of bandwidths that
#' will be applied to each function, or a matrix with individual bandwidths applied to
#' learning functions, one row per function. Alternatively, if \code{kNN} is specified,
#' the bandwidth matrix is built automatically using the local nearest neighbors bandwidths.
#' @param kNN       sequence of nearest neighbors (as a fraction of total number of neighbors, i.e.
#' a number between 0 and 1) to be considered for local bandwidth
#' selection. If \code{kNN} is provided, matrix \code{H} is build automatically by
#' evaluating sample quantiles of the distances in rows of the matrix of distances \code{D}.
#' @param nCV       if provided, only the first \code{nCV} leave-one-out cross-validated
#' trials is considered in the resulting mean squared error. By default, all functions are
#' included in cross-validation.
#' @param kernel    kernel function. Possible choices are \code{"uniform"}, \code{"triangular"},
#' \code{"Epanechnikov"}, \code{"biweight"}, \code{"triweight"}, and \code{"Gaussian"}.
#' @param quantile.type if \code{kNN} is provided, parameter that enters function \code{\link{quantile}}
#' that specifies the type of the sample quantile taken
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{CV} A matrix that corresponds to \code{H} that contains all the resulting
#' mean squared errors of the cross-validation procedure.
#' \item \code{CVB} A matrix that corresponds to \code{H} that contains all the resulting
#' variance estimates (mean squared errors - bias squared) of the cross-validation procedure.
#' \item \code{kernel} Kernel used.
#' }
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#'	@references Ferraty, F., and Nagy, S. (2021).
#'	Scalar-on-function local linear regression and beyond.
#'	\emph{Biometrika}, to appear.
#'
#' @examples
#' dat = generate(100,200,4)
#' LEARN = dat$LEARN
#' PRED = dat$PRED
#' Responses = dat$Responses
#' BASIS = dat$BASIS
#'
#' C = project(LEARN,BASIS)
#' D = L2metric(LEARN,LEARN)
#' llsm_cv(C, D, Responses, kNN = seq(.1,.9,length=11), kernel = "Gauss")
#'

llsm_cv <- function(C,D,Y,H=NULL,kNN=NULL,nCV=NULL, kernel=c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian"), quantile.type=7) {

  krn = c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian")
  kernel = match.arg(kernel)
  kernI = which(kernel==krn)
  if(is.null(H) & is.null(kNN)) stop("Neither candidate bandwidths nor nearest neighbour proportions specified.")
  if((!is.null(H)) & (!is.null(kNN))) stop("Specify either bandwidths or nearest neighbour proportions, not both.")
  # construct full matrix of candidate bandwidths, different for each function
  if(is.null(H)) Heval = length(kNN)
  if(is.null(kNN)){ if(is.null(dim(H))) Heval = length(H) else Heval = ncol(H) }
  Hf = matrix(nrow=nrow(D),ncol=Heval) # each function has now its own sequence of candidate bw
  for(i in 1:nrow(D)){
    if(is.null(H)){
      Hf[i,] = quantile(D[i,-i],kNN,type=quantile.type)
      Hf[i,kNN==1] = Inf
    }
    if(is.null(kNN)){
      if(is.null(dim(H))) Hf[i,] = H else Hf[i,] = H[i,]
    }
  }
  Hind = matrix(0,nrow=nrow(D), ncol=Heval);
  # Hind = (Hf>max(apply(D,1,function(x) quantile(x,10/nrow(D))))) # indicator for which H to perform also derivative cross-validation
  Hf[!is.finite(Hf)] = -1  # -1 is the code for usual linear regression
  if(is.null(C)){
    iC = matrix(1,nrow=nrow(D))
    iCnew = 1
  } else {
    iC<-cbind(1,C)
  }
  if(is.null(nCV)) nCV = nrow(D)
  if(nCV>nrow(D)) nCV = nrow(D)
  out <- .C("llsm_cv",
            C = as.double(iC), D = as.double(D), Y = as.double(Y),
            n = as.integer(nrow(D)), J = as.integer(ncol(iC)),
            Hf = as.double(Hf), Hind = as.integer(Hind), nH = as.integer(Heval),
            CV = rep(0,Heval*ncol(iC)),
            CVB = rep(0,Heval*ncol(iC)),
            CVD = rep(0,Heval*ncol(iC)),
            nCV = as.integer(nCV),
            kernI = as.integer(kernI),
            res = as.double(rep(0,ncol(iC)*ncol(iC)))
            # Yt = as.double(rep(0,ncol(iC)*ncol(iC))),
            # Ci = as.double(rep(0,(nrow(D)-1)*ncol(iC))),
            # Yi = as.double(rep(0,nrow(D)-1)),
            # Di = as.double(rep(0,nrow(D)-1)),
            # Cinew = as.double(rep(0,ncol(iC)))
  )
  CVM = matrix(out$CV,nrow=Heval)
  CVM[is.na(CVM)] = Inf
  CVB = matrix(out$CVB,nrow=Heval)
  CVB[is.na(CVB)] = Inf
  CVDM = matrix(out$CVD,nrow=Heval)
  CVDM[is.na(CVDM)] = Inf
  # if(inf){  # if Inf is in H, perform usual linear regression too (limit case as h->Inf)
  #   # that is the same as to ignore the information about the distances (D)
  #   # for derivative-based cross-validation, take H very large (matters only in this cross-validation procedure)
  #   R = llsm_cv(C,D,Y,H=Inf,nCV = nCV)
  #   CVM = rbind(CVM,R$CV)
  #   CVDM = rbind(CVDM,R$CVD)
  #   H = c(H,Inf)
  # }
  Hf[Hf==-1] = Inf # back to the usual coding for Inf
  if(!is.null(H)){ # H cross-validation
    # rownames(CVM) = rownames(CVB) = rownames(CVDM) = paste("H=",round(H,2),sep="")
    colnames(CVM) = colnames(CVB) = colnames(CVDM) = paste("J=",0:ncol(C),sep="")
    # h = H[apply(CVM,2,which.min)]
    # hD = H[apply(CVDM,2,which.min)]
    return(list(CV = CVM, CVB = CVB, kernel = kernel))
  }
  if(is.null(H)){ # kNN cross-validation
    rownames(CVM) = rownames(CVB) = rownames(CVDM) = paste("kNN=",round(kNN,2),sep="")
    colnames(CVM) = colnames(CVB) = colnames(CVDM) = paste("J=",0:ncol(C),sep="")
    # h = H[apply(CVM,2,which.min)]
    # hD = H[apply(CVDM,2,which.min)]
    return(list(CV=CVM, CVB = CVB, kernel = kernel))
  }
}

llsm_cv_single <- function(C,D,Y,H=NULL,kNN=NULL,nCV=NULL, kernel=c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian"), quantile.type=7) {
  ## input:
  # C:      basis coefficients of the X functions (matrix n*J)
  # D:      matrix of distances of X from X functions (matrix n*n)
  # Y:      response vector (vector n)
  # H:      sequence of bandwidths to cross-validate for
  # in CV for derivative, only such H are taken so that at least 10 functions are in
  # the h-neighborhood of any function x / otherwise we get very unstable, unreliable
  # results
  # kNN:    proportion of k-nearest neighbors to be used in bandwidth definition (in the interval (0,1))
  # quantile.type: type option in quantile function (see ?quantile)

  krn = c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian")
  kernel = match.arg(kernel)
  kernI = which(kernel==krn)
  if(is.null(H) & is.null(kNN)) stop("Neither candidate bandwidths nor nearest neighbour proportions specified.")
  if((!is.null(H)) & (!is.null(kNN))) stop("Specify either bandwidths or nearest neighbour proportions, not both.")
  # construct full matrix of candidate bandwidths, different for each function
  if(is.null(H)) Heval = length(kNN)
  if(is.null(kNN)){ if(is.null(dim(H))) Heval = length(H) else Heval = ncol(H) }
  Hf = matrix(nrow=nrow(D),ncol=Heval) # each function has now its own sequence of candidate bw
  for(i in 1:nrow(D)){
    if(is.null(H)){
      Hf[i,] = quantile(D[i,-i],kNN,type=quantile.type)
      Hf[i,kNN==1] = Inf
    }
    if(is.null(kNN)){
      if(is.null(dim(H))) Hf[i,] = H else Hf[i,] = H[i,]
    }
  }
  Hind = matrix(0,nrow=nrow(D), ncol=Heval);
  # Hind = (Hf>max(apply(D,1,function(x) quantile(x,10/nrow(D))))) # indicator for which H to perform also derivative cross-validation
  Hf[!is.finite(Hf)] = -1  # -1 is the code for usual linear regression
  if(is.null(C)){
    iC = matrix(1,nrow=nrow(D))
    iCnew = 1
  } else {
    iC<-cbind(1,C)
  }
  if(is.null(nCV)) nCV = nrow(D)
  if(nCV>nrow(D)) nCV = nrow(D)
  out <- .C("llsm_cv_single", 
            C = as.double(iC), D = as.double(D), Y = as.double(Y),
            n = as.integer(nrow(D)), J = as.integer(ncol(iC)),
            Hf = as.double(Hf), # Hind = as.integer(Hind),
            nH = as.integer(Heval),
            CV = rep(0,Heval*ncol(iC)),
            CVB = rep(0,Heval*ncol(iC)),
            # CVD = rep(0,Heval*ncol(iC)),
            nCV = as.integer(nCV),
            kernI = as.integer(kernI),
            res = as.double(rep(0,ncol(iC)*ncol(iC)))
            # Yt = as.double(rep(0,ncol(iC)*ncol(iC))),
            # Ci = as.double(rep(0,(nrow(D)-1)*ncol(iC))),
            # Yi = as.double(rep(0,nrow(D)-1)),
            # Di = as.double(rep(0,nrow(D)-1)),
            # Cinew = as.double(rep(0,ncol(iC)))
  )
  CVM = matrix(out$CV,nrow=Heval)
  CVM[is.na(CVM)] = Inf
  CVB = matrix(out$CVB,nrow=Heval)
  CVB[is.na(CVB)] = Inf
  # if(inf){  # if Inf is in H, perform usual linear regression too (limit case as h->Inf)
  #   # that is the same as to ignore the information about the distances (D)
  #   # for derivative-based cross-validation, take H very large (matters only in this cross-validation procedure)
  #   R = llsm_cv(C,D,Y,H=Inf,nCV = nCV)
  #   CVM = rbind(CVM,R$CV)
  #   CVDM = rbind(CVDM,R$CVD)
  #   H = c(H,Inf)
  # }
  Hf[Hf==-1] = Inf # back to the usual coding for Inf
  if(!is.null(H)){ # H cross-validation
    # rownames(CVM) = rownames(CVB) = rownames(CVDM) = paste("H=",round(H,2),sep="")
    colnames(CVM) = colnames(CVB) = paste("J=",0:ncol(C),sep="")
    # h = H[apply(CVM,2,which.min)]
    # hD = H[apply(CVDM,2,which.min)]
    return(list(CV = CVM, CVB = CVB, kernel = kernel))
  }
  if(is.null(H)){ # kNN cross-validation
    rownames(CVM) = rownames(CVB) = paste("kNN=",round(kNN,2),sep="")
    colnames(CVM) = colnames(CVB) = paste("J=",0:ncol(C),sep="")
    # h = H[apply(CVM,2,which.min)]
    # hD = H[apply(CVDM,2,which.min)]
    return(list(CV=CVM, CVB = CVB, kernel = kernel))
  }
}


lqsm_cv_single <- function(C,D,Y,H=NULL,kNN=NULL,nCV=NULL, kernel=c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian"), quantile.type=7) {
  ## input:
  # C:      basis coefficients of the X functions (matrix n*1)
  # D:      matrix of distances of X from X functions (matrix n*n)
  # Y:      response vector (vector n)
  # H:      sequence of bandwidths to cross-validate for
  # in CV for derivative, only such H are taken so that at least 10 functions are in
  # the h-neighborhood of any function x / otherwise we get very unstable, unreliable
  # results
  # kNN:    proportion of k-nearest neighbors to be used in bandwidth definition (in the interval (0,1))
  # quantile.type: type option in quantile function (see ?quantile)

  krn = c("uniform","triangular","Epanechnikov","biweight","triweight","Gaussian")
  kernel = match.arg(kernel)
  kernI = which(kernel==krn)
  if(is.null(H) & is.null(kNN)) stop("Neither candidate bandwidths nor nearest neighbour proportions specified.")
  if((!is.null(H)) & (!is.null(kNN))) stop("Specify either bandwidths or nearest neighbour proportions, not both.")
  # construct full matrix of candidate bandwidths, different for each function
  if(is.null(H)) Heval = length(kNN)
  if(is.null(kNN)){ if(is.null(dim(H))) Heval = length(H) else Heval = ncol(H) }
  Hf = matrix(nrow=nrow(D),ncol=Heval) # each function has now its own sequence of candidate bw
  for(i in 1:nrow(D)){
    if(is.null(H)){
      Hf[i,] = quantile(D[i,-i],kNN,type=quantile.type)
      Hf[i,kNN==1] = Inf
    }
    if(is.null(kNN)){
      if(is.null(dim(H))) Hf[i,] = H else Hf[i,] = H[i,]
    }
  }
  Hind = matrix(0,nrow=nrow(D), ncol=Heval);
  # Hind = (Hf>max(apply(D,1,function(x) quantile(x,10/nrow(D))))) # indicator for which H to perform also derivative cross-validation
  Hf[!is.finite(Hf)] = -1  # -1 is the code for usual linear regression
  if(is.null(C)){
    iC = matrix(1,nrow=nrow(D))
    iCnew = 1
  } else {
    iC<-cbind(1,C)
  }
  if(is.null(nCV)) nCV = nrow(D)
  if(nCV>nrow(D)) nCV = nrow(D)
  out <- .C("lqsm_cv_single",
            C = as.double(iC), D = as.double(D), Y = as.double(Y),
            n = as.integer(nrow(D)), J = as.integer(ncol(iC)),
            Hf = as.double(Hf), # Hind = as.integer(Hind),
            nH = as.integer(Heval),
            CV = rep(0,Heval*ncol(iC)),
            CVB = rep(0,Heval*ncol(iC)),
            # CVD = rep(0,Heval*ncol(iC)),
            nCV = as.integer(nCV),
            kernI = as.integer(kernI),
            res = as.double(rep(0,ncol(iC)*ncol(iC)))
            # Yt = as.double(rep(0,ncol(iC)*ncol(iC))),
            # Ci = as.double(rep(0,(nrow(D)-1)*ncol(iC))),
            # Yi = as.double(rep(0,nrow(D)-1)),
            # Di = as.double(rep(0,nrow(D)-1)),
            # Cinew = as.double(rep(0,ncol(iC)))
  )
  CVM = matrix(out$CV,nrow=Heval)
  CVM[is.na(CVM)] = Inf
  CVB = matrix(out$CVB,nrow=Heval)
  CVB[is.na(CVB)] = Inf
  # if(inf){  # if Inf is in H, perform usual linear regression too (limit case as h->Inf)
  #   # that is the same as to ignore the information about the distances (D)
  #   # for derivative-based cross-validation, take H very large (matters only in this cross-validation procedure)
  #   R = llsm_cv(C,D,Y,H=Inf,nCV = nCV)
  #   CVM = rbind(CVM,R$CV)
  #   CVDM = rbind(CVDM,R$CVD)
  #   H = c(H,Inf)
  # }
  Hf[Hf==-1] = Inf # back to the usual coding for Inf
  if(!is.null(H)){ # H cross-validation
    # rownames(CVM) = rownames(CVB) = rownames(CVDM) = paste("H=",round(H,2),sep="")
    colnames(CVM) = colnames(CVB) = paste("J=",0:ncol(C),sep="")
    # h = H[apply(CVM,2,which.min)]
    # hD = H[apply(CVDM,2,which.min)]
    return(list(CV = CVM, CVB = CVB, kernel = kernel))
  }
  if(is.null(H)){ # kNN cross-validation
    rownames(CVM) = rownames(CVB) = paste("kNN=",round(kNN,2),sep="")
    colnames(CVM) = colnames(CVB) = paste("J=",0:ncol(C),sep="")
    # h = H[apply(CVM,2,which.min)]
    # hD = H[apply(CVDM,2,which.min)]
    return(list(CV=CVM, CVB = CVB, kernel = kernel))
  }
}

#' Basis projection
#'
#' Project a set of functions onto a basis.
#'
#' @param X         a vector of a single function, or a matrix with one function per row
#' @param basis     a matrix of basis functions, one function per row
#'
#' @return A matrix of coefficients of the functions in \code{X} with respect to the
#' basis in \code{basis}.
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#' @examples
#' gridsize = 101                           # grid size
#' Grid = seq(0, 1, length = gridsize)      # grid of measurements
#' n = 150
#' Jmodel = 4                               # number of basis elements
#' BASIS = phif(Jmodel, Grid)               # first Jmodel Fourier basis elements
#' UNIF = matrix(runif(sample.size * Jmodel, min = -1, max = 1), nrow = sample.size, ncol = Jmodel)
#'
#' PREDICTORS = UNIF %*% BASIS[1:Jmodel, ]
#' PREDICTORS2 = expand(UNIF,BASIS[1:Jmodel, ])
#' all.equal(PREDICTORS,PREDICTORS2)        # expansion of coefficients into functions
#'
#' COEF = project(PREDICTORS,BASIS)         # projection back to basis functions
#' max(abs(COEF-UNIF))

project = function(X,basis){
  # projects functions in rows of X onto the basis given by the rows of basis
  if(is.vector(X)) X = matrix(X,nrow=1)
  if(is.vector(basis)) basis = matrix(basis,nrow=1)
  n = nrow(X); d = ncol(X); J = nrow(basis);
  if(ncol(basis)!=d) stop("project: dimension of the functions does not match")
  C = matrix(nrow=n, ncol=J) # matrix of coefficients
  for(i in 1:n) C[i,] = apply(basis,1,function(x) mean(x*X[i,]))
  return(C)
}

#' Basis expansion
#'
#' Expand a set of coefficients into functions in a basis.
#'
#' @param C         a vector of coefficients, or a matrix with one vector per row
#' @param basis     a matrix of basis functions, one function per row
#'
#' @return A matrix of functional values that corresponds to the coefficients given
#' in the rows of \code{C} with respect to the basis in \code{basis}.
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#' @examples
#' gridsize = 101                           # grid size
#' Grid = seq(0, 1, length = gridsize)      # grid of measurements
#' n = 150
#' Jmodel = 4                               # number of basis elements
#' BASIS = phif(Jmodel, Grid)               # first Jmodel Fourier basis elements
#' UNIF = matrix(runif(n * Jmodel, min = -1, max = 1), nrow = n, ncol = Jmodel)
#'
#' PREDICTORS = UNIF %*% BASIS[1:Jmodel, ]
#' PREDICTORS2 = expand(UNIF,BASIS[1:Jmodel, ])
#' all.equal(PREDICTORS,PREDICTORS2)        # expansion of coefficients into functions
#'

expand = function(C,basis,na.to.zero=TRUE){
  # expands the coefficient from the rows of C into the basis given by the rows of basis
  if(is.vector(C)) C = matrix(C,nrow=1)
  if(is.vector(basis)) basis = matrix(basis,nrow=1)
  n = nrow(C); J = ncol(C);
  # if(ncol(basis)!=d) stop("expand: dimension of basis functions is not d")
  if(J!=nrow(basis)) stop("expand: dimension of coefficients does not match the number of basis functions")
  if(na.to.zero) C[is.na(C)] = 0
  X = C%*%basis
  return(X)
}

#' Fourier basis functions
#'
#' Creates a matrix of the first few functions from the Fourier basis.
#'
#' @param J Number of basis functions.
#' @param t Vector of the common points in the domain where the basis functions are evaluated.
#'
#' @return A matrix with \code{J} rows. In each row there are functional values of a basis
#' function that correspond to the points in the domain from \code{t}.
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#' @examples
#' gridsize = 101                           # grid size
#' J = 4                                    # number of basis elements
#' BASIS = phif(Jmodel, Grid)               # first J Fourier basis elements

phif = function(J,t){
  # trigonometric basis of J functions
  # output is matrix of dimensions (J*length(t))
  b = (J-1)/2+1
  B = matrix(nrow=2*b+1,ncol=length(t))
  B[1,] = 1
  for(j in 1:b){
    B[2+2*(j-1),] = sqrt(2)*sin(2*pi*j*t)
    B[3+2*(j-1),] = sqrt(2)*cos(2*pi*j*t)
  }
  return(B[1:J,])
}

# #' B-spline basis functions
# #'
# #' Creates a matrix of the first few functions from the B-spline basis without knots.
# #'
# #' @param J Number of basis functions.
# #' @param t Vector of the common points in the domain where the basis functions are evaluated.
# #'
# #' @return A matrix with \code{J} rows. In each row there are functional values of a basis
# #' function that correspond to the points in the domain from \code{t}.
# #'
# #'
# #' @examples
# #' gridsize = 101                           # grid size
# #' J = 4                                    # number of basis elements
# #' BASIS = phib(Jmodel, Grid)               # first J basis elements
#
# phib = function(J,t) t(bs(t,df=J,intercept=TRUE))
# # B-spline basis of J functions, without knots, matrix of dimensions (J*length(t))

#' \code{L2}-metric
#'
#' Fast computation of the \code{L2}-metric between two sets of functional data.
#'
#' @param A A matrix of functional values, one function per row, \code{n} rows.
#' @param B A matrix of functional values, one function per row, \code{m} rows.
#' @param correct Logical; should the boundary points be computed only by a half of their
#' value due to boundary effects?
#'
#' @return A matrix of size \code{n}-times-\code{m} that contains approximations to the
#' \code{L2}-distances between function that correspond to the rows of \code{A} and \code{B}.
#' The functions must be evaluated at an equidistant grid. It may be necessary to multiply
#' by a constant that corresponds to the grid size.
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#' @examples
#' gridsize = 101                           # grid size
#' J = 4                                    # number of basis elements
#' BASIS = phif(Jmodel, Grid)
#' L2metric(BASIS, BASIS)
#' L2metric(BASIS, BASIS, correct = FALSE)
#' sqrt(crossprod(BASIS[1,]-BASIS[2,]))
#'

L2metric = function(A,B,correct=TRUE){
  # computes fast approximation of L2 distance between fctions A and B
  if(correct){
    M = .Fortran("metrl2",
                 as.numeric(A),
                 as.numeric(B),
                 as.integer(m<-nrow(A)),
                 as.integer(n<-nrow(B)),
                 as.integer(d<-length(as.numeric(A))/nrow(A)),
                 m = as.numeric(rep(-1,m*n)))$m
  } else {
    M = .Fortran("metrl2b",
                 as.numeric(A),
                 as.numeric(B),
                 as.integer(m<-nrow(A)),
                 as.integer(n<-nrow(B)),
                 as.integer(d<-length(as.numeric(A))/nrow(A)),
                 m = as.numeric(rep(-1,m*n)))$m
  }
  return(M = matrix(M,nrow=m))
}
