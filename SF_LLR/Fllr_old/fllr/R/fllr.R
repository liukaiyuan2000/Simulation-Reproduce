################################################################################
# Functional local linear regression (fllr) including data-driven basis
# (functional predictor's FPCA) and an automatic choice of:
#  * dimension J of the approximating subspace
#  * bandwidth for the regression operator
#  * bandwidth for the functional derivative
################################################################################

#' Functional local linear regression
#'
#' Functional local linear regression (fllr) with automatic bandwidth selection
#' including data-driven basis (functional predictor's FPCA) and a data driven
#' selection of the dimension of the approximating subspace.
#'
#' @param Responses vector of scalar responses
#' @param LEARN     matrix of regressors, functional values evaluated at common
#' domain, one function per row
#' @param PRED      matrix of functions at which new values are predicted, one
#' function per row
#' @param nboot     number of wild bootstrap replicates used for finding the
#' bandwidth for functional derivatives
#' @param kNN.grid.length  grid length for the bandwidth search
#' @param percent   if \code{BASIS} is not given, percentage of explained sample
#' variance that is taken into account when building the basis expansion in
#' terms of eigenfunctions of the empirical covariance operator
#' @param boot.seed random seed for the bootstrap procedure
#' @param BASIS     if basis functions are not estimated using the empirical
#' covariance operator, matrix of basis functions. One function per row.
#' @param method    bandwidth selection method for the regression operator.
#' Possible values are \code{cv} for leave-one-out cross-validation, and
#' \code{aicc} for the improved Akaike information criterion
#' @param derivative logical; should also the bandwidth for the functional
#' derivative be searched for? If set to \code{FALSE}, the bandwidth for the
#' regression operator is used also for the functional derivatives.
#' @param J.est      logical; should it be cross-validated also with respect to
#' the dimension of the basis? If set to \code{FALSE}, the number of functions
#' in \code{BASIS} is taken to be the optimal dimension. Ignored if
#' \code{BASIS = NULL}.
#' @param fullgrid   logical; should the full grid [0,1] of k-NN be used
#' directly for cross-validation? If \code{FALSE}, the grid is built iteratively
#' from small to larger \code{k}. By default, \code{fullgrid=TRUE} if the number
#' of functions in \code{LEARN} is at most 250.
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{Estimated.responses} A vector of responses estimated by the local
#' linear method that correspond to the functions in the rows of \code{LEARN}.
#' \item \code{ESTIMATED.DERIV} A matrix of the estimated (Riesz representations
#' of) functional derivatives of the regression operator, each row corresponding
#' to a row of \code{LEARN}.
#' \item \code{Predicted.responses} If matrices \code{LEARN} and \code{PRED}
#' differ, a vector of predicted responses, each element of the vector
#' corresponding to a row of \code{PRED}.
#' \item \code{PREDICTED.DERIV} If matrices \code{LEARN} and \code{PRED} differ,
#' a matrix of predicted (Riesz representations of) functional derivatives,
#' each row corresponding to a row of \code{PRED}.
#' \item \code{min.crit.reg} Minimum criterion value (mean squared error for
#' \code{method = cv} or AICc for \code{method = aicc}) obtained in the search
#' for the optimal bandwidth for the regression operator.
#' \item \code{knn.opt} List of optimal bandwidths for the regression operator
#' \code{reg} and the functional derivative \code{deriv}.
#' \item \code{dim.opt} Optimal basis dimension used for both the regression
#' operator and the functional derivative.
#' }
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#' @references Ferraty, F., and Nagy, S. (2021).
#'	Scalar-on-function local linear regression and beyond.
#'	\emph{Biometrika}, to appear.
#'
#' @examples
#' dat = generate(100,200,4)
#' LEARN = dat$LEARN
#' PRED = dat$PRED
#' Responses = dat$Responses
#'
#' res = fllr(Responses, LEARN, PRED)
#'
#' plot(res$Predicted.responses~dat$oracle.Responses)
#'

fllr = function(Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30,
                percent = 0.95, boot.seed = NULL, BASIS = NULL, method = "cv",
                derivative = TRUE, J.est = TRUE, fullgrid=NULL)
################################################################################
# The "fllr" routine estimates the Functional Local Linear Regression as well as
# its functional derivative. "fllr" computes a data-driven basis (functional
# predictor's FPCA) and selects the dimension of the approximating subspace.
# Two specific bandwidths are automatically chosen: one for estimating the
# regression operator and one for the functional derivative.
################################################################################
{
  if(is.null(fullgrid)) fullgrid = (nrow(LEARN)<=250)
  # Check whether there are additional observations for predicting
  # corresponding responses
  outsample = T
  if(identical(LEARN,PRED)) outsample = F
  ## if LEARN is the same matrix as PRED, outsample = F and PRED = LEARN
  ##
  # Estimate the second order moment operator
  gridsize = ncol(LEARN)
  if(is.null(BASIS)){
    if(outsample){
      sample.size = nrow(LEARN) + nrow(PRED)
      # COV = crossprod(rbind(LEARN, PRED)) / (sample.size * (gridsize - 1))
      COV = var(rbind(LEARN,PRED))/(gridsize-1) # added in R1
    }else{
      sample.size = nrow(LEARN)
      # COV = crossprod(LEARN) / (sample.size * (gridsize - 1))
      COV = var(LEARN)/(gridsize-1) # added in R1
    }
  }
  # For a possibly slightly faster version
  # Install the r package "parallelDist" to use "parDist" instead of "dist"
  # to compute distances ||X_i - X_j||^2 between functional predictors
  # require(parallelDist)
  # Dependency on package "parallelDist" however causes problems in R version 4
  DIST = as.matrix(dist(LEARN, method = "euclidean")) / sqrt(gridsize - 1)
  DIST.SORTED.ROW.BY.ROW = t(apply(DIST, 1, sort))
  if(is.null(BASIS)){
    # Compute estimated first "dim.max" eigenfunctions
    res.eig = eigen(COV, sym = T)
    # Drop first eigenvalue and derive the approximating subspace expressing at
    # least "percent" of the remaining total variance
    cum.part.of.variance = cumsum(res.eig$values[-1]) / sum(res.eig$values[-1])
    dim.min = 1
    dim.max = min(sum(cum.part.of.variance < percent) + 2,gridsize)
    BASIS = t(res.eig$vectors[, 1:dim.max]) * sqrt(gridsize - 1)
  } else {
    # if J is estimated, dimensions are searched from one
    if(J.est) dim.min = 1 else dim.min = nrow(BASIS)
    dim.max = nrow(BASIS)
  }
  nlearn = nrow(LEARN)
  INNER.LEARN = tcrossprod(LEARN, BASIS) / (gridsize - 1)
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
    ############################################################################
    # estimate regression operator with bandwidth maximizing cross-validation
    # criterion
    ############################################################################
    CRIT.REG = matrix(Inf, actual.kNN.grid.length, dim.max)
    if(method=="cv"){
      H = matrix(nrow=nrow(INNER.LEARN),ncol=actual.kNN.grid.length)
      for(knni in 1:actual.kNN.grid.length)
        H[,knni] = DIST.SORTED.ROW.BY.ROW[, Knn.grid[knni] + 1]
      for(dim in dim.min:dim.max){
        CRIT.REG[,dim] = suppressWarnings(
          llsm_cv_single(INNER.LEARN[,1:dim,drop=FALSE],
                         DIST,Responses,H=H,kernel="Epanechnikov")$CV[,dim+1])
        # warnings indicate that the Cholesky decomposition failed, in case
        # when kNN is smaller than the dimension of the model. Then, the CV
        # criteria are are replaced by Inf
        CRIT.REG[Knn.grid<dim+2,dim] = Inf
        # numerical fix
        # in the case knn<dim, the weight matrix had to be singular for kernels
        # supported in [0,1]
      }
      dim.opt = (1:dim.max)[which.min(apply(CRIT.REG,2,min))]
      # plot(apply(CRIT.REG,2,min),type="b")
    }
    if(method=="aicc"){
      for(dim in (dim.min:dim.max)){
        # compute inner products between functional predictors and basis
        # functions
        INNER.LEARN2 = INNER.LEARN[,1:dim,drop=FALSE]
        Estimated.responses = 0
        count = 0
        for(knn in Knn.grid){
          count = count + 1
          trHatMat = 0
          for(ii in 1:nlearn){
            PHI = cbind(1, t(t(INNER.LEARN2) - INNER.LEARN2[ii, ]))
            bwd = DIST.SORTED.ROW.BY.ROW[ii, knn + 1]
            U = DIST[ii,] / bwd
            Kernel = (1 - U^2) * (U < 1)
            PHIK = PHI * Kernel
            TEMP = crossprod(PHIK, PHI)
            TEMP.INV = tryCatch(solve(TEMP),
                         error=function(x) return(matrix(Inf,nrow=nrow(TEMP),
                                                         ncol=ncol(TEMP))))
            Temp = crossprod(PHIK, Responses)
            Estimated.responses[ii] = crossprod(TEMP.INV[1, ], Temp)
            # Compute the trace of the hat matrix for the AICc criterion
            trHatMat = trHatMat + crossprod(TEMP.INV[1, ], t(PHIK)[, ii])
          }
          # In case the matrix TEMP was not invertible (we obtained a matrix
          # of Inf's in TEMP.INV) we set CRIT.REG = Inf
          # Compute criterion for selecting the optimal bandwidth
          Estimated.errors = Estimated.responses - Responses
          CRIT.REG[count, dim] = mean(Estimated.errors^2)
          CRIT.REG[count, dim] = log(CRIT.REG[count, dim]) + 1 +
            (2 * trHatMat + 2)/(nlearn - trHatMat - 2)
          if(is.na(CRIT.REG[count, dim])) CRIT.REG[count,dim] = Inf
          if(knn<dim+2) CRIT.REG[count,dim] = Inf
          # numerical fix, take at least dim+2 nearest neighbors to obtain
          # numerically stable results
        }
        dim.opt = dim
        if(dim >= 2){
          if(min(CRIT.REG[, dim - 1]) < min(CRIT.REG[, dim])){
            dim.opt = dim - 1
            break
          }
        }
      }
    }
    if(rank.min != order(CRIT.REG[, dim.opt])[1]) break
  }
  # Determine corresponding criterion value for the regression
  min.crit.reg = min(CRIT.REG[, dim.opt])
  # Determine optimal number of k-nearest neighbours and dimension J for the
  # regression
  knn.reg = (Knn.grid)[order(CRIT.REG[, dim.opt])[1]]
  # Compute inner products between functional predictors and basis functions
  B = matrix(BASIS[1:dim.opt,],nrow=dim.opt)
  INNER.LEARN = matrix(INNER.LEARN[,1:dim.opt],ncol=dim.opt)
  if(method=="cv") Estimated.responses = llsm_leave(
    INNER.LEARN,INNER.LEARN,DIST,Responses,
    h = DIST.SORTED.ROW.BY.ROW[, knn.reg + 1], kernel="Epanechnikov")[,1]
  if(method=="aicc") Estimated.responses = llsm(
    INNER.LEARN,INNER.LEARN,DIST,Responses,
    h = DIST.SORTED.ROW.BY.ROW[, knn.reg + 1], kernel="Epanechnikov")[,1]
  ##############################################################################
  # estimate fuctional derivatives (with specific bandwidth choice)
  ##############################################################################
  if(derivative){
    # Step 1: compute estimated errors
    Estimated.errors = Responses - Estimated.responses
    # Set values for the wild bootstrap procedure
    gn = 0.5 * (1 + sqrt(5))
    gn.conj = 0.5 * (1 - sqrt(5))
    prob = (5 + sqrt(5)) / 10
    knn.offset = 0 # offset of the knn bandwidth for pilot estimator,
                   # increases in case of numerical instabilities
    stable = FALSE # indicator whether the search for the pilot estimator
                   # is numerically stable, if not, knn.offset increases
    while(!stable){
    # Initialization
    Crit.deriv.boot = 0
    COEF.DERIV.BOOT = matrix(0, nlearn, dim.opt)
    PHIK = list()
    TEMP.INV = list()
    # Compute invariant quantities
    for(ii in 1:nlearn){
      PHI = cbind(1, t(t(INNER.LEARN) - INNER.LEARN[ii, ]))
      bwd = DIST.SORTED.ROW.BY.ROW[ii, knn.reg + 1 + knn.offset]
      U = DIST[ii, 1:nlearn] / bwd
      Kernel = (1 - U^2) * (U < 1)
      Kernel[ii] = 0
      PHIK[[ii]] = PHI * Kernel
      TEMP = crossprod(PHIK[[ii]], PHI)
      TEMP.INV[[ii]] = tryCatch(solve(TEMP),
                        error=function(x) return(matrix(Inf,nrow=nrow(TEMP),
                                                            ncol=ncol(TEMP))))
    }
    # Repeat B=100 (by default) times the wild bootstrap procedure
    if(!is.null(boot.seed)) set.seed(boot.seed)
    for(bb in 1:nboot){
      # Step 2: build bootstrapped errors
      U  = runif(nlearn)
      V = (U >= prob) * gn + (U < prob) * gn.conj
      Boostrapped.errors = Estimated.errors * V
      # Step 3: drawn the bootstrap sample
      Boostrapped.responses = Estimated.responses + Boostrapped.errors
      # Step 4: compute bootstrapped functional derivatives
      for(ii in 1:nlearn){
        Temp = crossprod(PHIK[[ii]], Boostrapped.responses)
        # Estimate simultaneously the fuctional derivatives basis expansion
        COEF.DERIV.BOOT[ii, ] = COEF.DERIV.BOOT[ii, ] +
          TEMP.INV[[ii]][-1, ] %*% Temp
        for(iii in 1:dim.opt) if(is.na(COEF.DERIV.BOOT[ii,iii]))
          COEF.DERIV.BOOT[ii, iii] = Inf
      }
    }
    if(sum(is.infinite(COEF.DERIV.BOOT)|is.na(COEF.DERIV.BOOT))>0)
      knn.offset = knn.offset + 1 else stable = TRUE
    }
    # Compute the average of the boostrapped estimations for each functional
    # predictor in the learning sample
    COEF.DERIV.BOOT = COEF.DERIV.BOOT / nboot
    # Compute the bootstrap pilot
    PILOT.BOOT.DERIV = COEF.DERIV.BOOT %*% B
    # Determine the optimal bandwidth for estimating the functional derivatives
    rank.min = kNN.grid.length
    if(fullgrid) Knn.learn.ratio = c(0.02,1) else
      Knn.learn.ratio = c(0.02, 0.25, 0.5, 0.75, 1)
    # Build an adaptative grid of kNN (number of knearest neighbours used for
    # determining local bandwidths). If the minimum of the criterion is reached
    # for the last bandwidth then the grid "Knn.grid" is updated
    kNN.iterate = TRUE # indicator stopping the iteration in kNN CV procedure
    knn.max = dim.opt  # initiate the kNN CV procedure, for less than dim.opt
    # the weight matrix is singular (for kernels on [0,1])
    ratio = 1
    while(kNN.iterate){
      knn.min = min(max(round(nlearn * Knn.learn.ratio[ratio]), knn.max),
                    nlearn-1)
      knn.max = min(max(round(nlearn * Knn.learn.ratio[ratio + 1]),
                        knn.min + kNN.grid.length  - 1),nlearn-1)
      ratio = ratio + 1
      if(ratio == length(Knn.learn.ratio)-1) kNN.iterate= FALSE # stop iteration
      if(knn.max == nlearn-1) kNN.iterate= FALSE                # stop iteration
      #    for(ratio in 1:(length(Knn.learn.ratio) - 1)){
      #      knn.min = max(round(nlearn * Knn.learn.ratio[ratio]), dim.max + 4)
      #      knn.max = round(nlearn * Knn.learn.ratio[ratio + 1])
      knn.range = knn.max - knn.min
      if(knn.range <= kNN.grid.length){
        Knn.grid = knn.min:knn.max
      }else{
        step = ceiling(knn.range / kNN.grid.length)
        Knn.grid = seq(from = knn.min, to = knn.max, by = step)
      }
      actual.kNN.grid.length = length(Knn.grid)
      rank.min = actual.kNN.grid.length
      ##########################################################################
      # estimate kNN
      ##########################################################################
      Crit.deriv.boot = rep(NA,actual.kNN.grid.length)
      count = 0
      for(knn in Knn.grid){
        count = count + 1
        COEF.DERIV = suppressWarnings(
          llsm_leave(INNER.LEARN,INNER.LEARN,DIST,Responses,h=
                       DIST.SORTED.ROW.BY.ROW[, knn + 1],kernel="Epa")[,-1])
        ESTIMATED.DERIV = COEF.DERIV %*% B
        # Compute the criterion for estimating the fuctional derivative
        Crit.deriv.boot[count] = mean((PILOT.BOOT.DERIV - ESTIMATED.DERIV)^2)
      }
      Crit.deriv.boot[is.na(Crit.deriv.boot)] = Inf
      # plot(Crit.deriv.boot~Knn.grid,type="b")
      if(rank.min != order(Crit.deriv.boot)[1]) break
    }
    # Determine the criterion value
    crit.deriv.boot = min(Crit.deriv.boot)
    # Determine the optimal number of k-nearest neighbours for the functional
    # derivative
    knn.deriv = (Knn.grid)[which.min(Crit.deriv.boot)]
  } else { # if the badnwidth for the derivative is just that from regression
    knn.deriv = knn.reg
  }
  # Save optimal bandwidths in terms of k-nearest neighbours into outputs
  knn.opt = list(reg = knn.reg, deriv = knn.deriv)
  COEF.DERIV = llsm(INNER.LEARN,INNER.LEARN,DIST,Responses,
                    h=DIST.SORTED.ROW.BY.ROW[, knn.deriv + 1],kernel="Epa")[,-1]
  ESTIMATED.DERIV = COEF.DERIV %*% B
  # Compute predicted responses + predicted fuctional derivative if "outsample"
  # is true (i.e. one has additional observations for the functional predictor)
  if(outsample){
    DIST = as.matrix(dist(rbind(LEARN, PRED), method = "euclidean")) / sqrt(
      gridsize - 1)
    DIST.SORTED.PRED.BY.LEARN = t(apply(DIST[-(1:nlearn), 1:nlearn], 1, sort))
    # Compute inner products between functional predictors and basis functions
    INNER.PRED = tcrossprod(PRED, B) / (gridsize - 1)
    nout = nrow(PRED)
    Predicted.responses = llsm(INNER.LEARN,INNER.PRED,DIST[-(1:nlearn),
                  1:nlearn],Responses,h=DIST.SORTED.PRED.BY.LEARN[,knn.reg + 1],
                  kernel="Epa")[,1]
    COEF.PRED.DERIV = llsm(INNER.LEARN,INNER.PRED,DIST[-(1:nlearn), 1:nlearn],
                           Responses,h=DIST.SORTED.PRED.BY.LEARN[,knn.deriv+ 1],
                           kernel="Epa")[,-1]
    PREDICTED.DERIV = COEF.PRED.DERIV %*% B
    return(list(Estimated.responses = Estimated.responses,
                ESTIMATED.DERIV = ESTIMATED.DERIV,
                Predicted.responses = Predicted.responses,
                PREDICTED.DERIV = PREDICTED.DERIV,
                min.crit.reg = min.crit.reg,
                knn.opt = knn.opt,
                dim.opt = dim.opt))
  }else{
    return(list(Estimated.responses = Estimated.responses,
                ESTIMATED.DERIV = ESTIMATED.DERIV,
                min.crit.reg = min.crit.reg,
                knn.opt = knn.opt,
                dim.opt = dim.opt))
  }
}
