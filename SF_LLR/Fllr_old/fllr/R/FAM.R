#' Functional additive model estimation
#'
#' Estimates the scalar-on-function model using the additive modeling approach
#' from Mueller and Yao (2008) for the regression operator, and Mueller and Yao
#' (2010) for the functional derivatives.
#'
#' @param Responses vector of scalar responses
#' @param LEARN     matrix of regressors, functional values evaluated at common
#' domain, one function per row
#' @param PRED      matrix of functions at which new values are predicted, one
#' function per row
#' @param method  method used for the selection of the number of eigenvectors,
#' possible values are \code{FVE} for the fraction-of-variance-explained in the
#' fitted covariance function, or \code{BIC} for the conditional pseudo-Gaussian
#' BIC selection
#' @param FVEthreshold  a numerical constant between 0 and 1, giving the
#' threshold for method \code{FVE}. Ignored with \code{method=BIC}.
#' @param BASIS     if basis functions are not estimated using the empirical
#' covariance operator, matrix of basis functions. One function per row.
#' @param knn.ind   logical; should the bandwidth be chosen using the k-NN
#' cross-validation? If \code{FALSE}, bandwidth is cross-validated using the
#' procedures from package \code{locpol}. The choice \code{TRUE} is compatible
#' with the bandwidth selection procedure in \code{fllr}.
#' @param fullgrid  logical; should the full grid [0,1] of k-NN be used
#' directly for cross-validation? If \code{FALSE}, the grid is built iteratively
#' from small to larger \code{k}. By default, \code{fullgrid=FALSE}.
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{Predicted.responses} A vector of predicted responses, each
#' element of the vector corresponding to a row of \code{PRED}.
#' \item \code{PREDICTED.DERIV} A matrix of predicted (Riesz representations
#' of) functional derivatives, each row corresponding to a row of \code{PRED}.
#' \item \code{PREDICTED.DERIV.LIN} Only if \code{knn.ind=TRUE}. A matrix of
#' predicted (Riesz representations of) functional derivatives, each row
#' corresponding to a row of \code{PRED}. Instead of in \code{PREDICTED.DERIV}
#' where the local quadratic estimator is employed, here the local linear
#' estimator is used. Only for compatibility purposes.
#' \item \code{dim} Dimension of the eigenbasis used for both the regression
#' operator and the functional derivative. Chosen as a given percentage of total
#' sample variance explained by the first eigenfunctions of the empirical
#' covariance operator (if \code{method}=\code{FVE}), or chosen using
#' the conditional pseudo-Gaussian BIC (if \code{method}=\code{BIC}).
#' }
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#' @references Mueller, H.-G., and Yao, F. (2008). Functional additive models.
#' \emph{Journal of the American Statistical Association} 103, 1534--1544.
#' @references Mueller, H.-G., and Yao, F. (2010). Additive modelling of
#' functional gradients. \emph{Biometrika} 97, 791--805.
#'
#' @examples
#' dat = generate(100,200,4)
#' LEARN = dat$LEARN
#' PRED = dat$PRED
#' Responses = dat$Responses
#'
#' res = FAM(Responses, LEARN, PRED)
#'
#' plot(res$Predicted.responses~dat$oracle.Responses)
#'

FAM = function(Responses, LEARN, PRED, method="FVE", FVEthreshold = .9,
              BASIS = NULL, knn.ind = TRUE, fullgrid=FALSE)
{
  outsample = T
  if(identical(LEARN,PRED)) outsample = F
  ## if LEARN is the same matrix as PRED, outsample = F
  # Estimate the second order moment operator
  gridsize = ncol(LEARN)
  if(is.null(BASIS)){
    if(outsample){
      # sample.size = nrow(LEARN) + nrow(PRED)
      # COV = crossprod(rbind(LEARN, PRED)) / (sample.size * (gridsize - 1))
      BASIS = eigenbasis(rbind(LEARN, PRED), method=method, FVEthreshold=
                           FVEthreshold)
    }else{
      # sample.size = nrow(LEARN)
      # COV = crossprod(LEARN) / (sample.size * (gridsize - 1))
      BASIS = eigenbasis(LEARN, method=method, FVEthreshold=FVEthreshold)
    }
    # Compute estimated first "dim.max" eigenfunctions
    # res.eig = eigen(COV, sym = T)
    # According to the guidelines from Mueller and Yao (2010, Section 3)
    # derive the approximating subspace expressing at least "percent"
    # of the total variance
    # cum.part.of.variance = cumsum(res.eig$values) / sum(res.eig$values)
    # dim = sum(cum.part.of.variance < percent) + 1
    # BASIS = t(res.eig$vectors[, 1:dim]) * sqrt(gridsize - 1)
    dim = ncol(BASIS)
    BASIS = t(BASIS) * sqrt(gridsize - 1)
  } else { dim = nrow(BASIS) }
  nlearn = nrow(LEARN)
  #
  B = BASIS
  INNER.LEARN = tcrossprod(LEARN, B) / (gridsize - 1)
  INNER.PRED = tcrossprod(PRED, B) / (gridsize - 1)
  nout = nrow(PRED)
  #
  if(!knn.ind) MY.raw = MuellerYao(INNER.LEARN,INNER.PRED,Responses)
  if(knn.ind)  MY.raw = MuellerYao_kNN(INNER.LEARN,INNER.PRED,Responses,
                                       fullgrid = fullgrid)
  if(!knn.ind){
    COEF.PRED.DERIV = MY.raw$coefficients
    PREDICTED.DERIV = COEF.PRED.DERIV %*% B
  }
  #
  if(knn.ind){
    COEF.PRED.DERIV = MY.raw$deriv.coefficients
    PREDICTED.DERIV = COEF.PRED.DERIV %*% B
    COEF.PRED.DERIV = MY.raw$deriv.coefficients.lin
    PREDICTED.DERIV.LIN = COEF.PRED.DERIV %*% B
  }
  #
  Predicted.responses = mean(Responses) + rowSums(MY.raw$reg.coefficients)
  if(knn.ind) return(list(Predicted.responses = Predicted.responses,
                          PREDICTED.DERIV = PREDICTED.DERIV,
                          PREDICTED.DERIV.LIN = PREDICTED.DERIV.LIN,
                          dim = dim)) else
              return(list(Predicted.responses = Predicted.responses,
                          PREDICTED.DERIV = PREDICTED.DERIV,
                          dim = dim))
}


MuellerYao_kNN = function(C, Cnew, Y, kNN.grid.length=30, fullgrid = FALSE){
  # additive functional derivative estimation (Mueller and Yao ,2010)
  ksi = C
  nlearn = nrow(C)
  nbasis = ncol(C)
  g = g1 = f = matrix(NA,nrow(Cnew),nbasis)
  Y = Y - mean(Y)     # centering Y
  for(k in 1:nbasis){
    # regression estimation
    X = ksi[,k]
    Xnew = Cnew[,k]
    DIST = abs(outer(X,X,"-"))
    X = matrix(X,ncol=1)
    Xnew = matrix(Xnew,ncol=1)
    DIST.SORTED.ROW.BY.ROW = t(apply(DIST, 1, sort))
    #
    if(nlearn<50) warning("With less than 50 learning curves the
                          cross-validation for the kNN method may be unstable.")
    if(fullgrid) Knn.learn.ratio = c(0.02,1) else
      Knn.learn.ratio = c(0.02, 0.25, 0.5, 0.75, 1)
    # Build an adaptative grid of kNN (number of knearest neighbours used for
    # determining local bandwidths). If the minimum of the criterion
    # ("CRIT.REG") is reached for the bandwidth corresponding to
    # "Knn.grid[kNN.grid.length]", then the grid "Knn.grid" is updated
    kNN.iterate = TRUE # indicator stopping the iteration in kNN CV procedure
    knn.max = 2        # initiate the kNN CV procedure
    ratio = 1
    while(kNN.iterate){
      knn.min = min(max(round(nlearn * Knn.learn.ratio[ratio]), knn.max),
                    nlearn-1)
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
      ##########################################################################
      # estimate regression operator with bandwidth maximizing the CV criterion
      ##########################################################################
      CRIT.REG = rep(NA,actual.kNN.grid.length)
      H = matrix(nrow=nlearn,ncol=actual.kNN.grid.length)
      for(knni in 1:actual.kNN.grid.length) H[,knni] = DIST.SORTED.ROW.BY.ROW[,
                                                             Knn.grid[knni] + 1]
      CRIT.REG = suppressWarnings(llsm_cv_single(X,DIST,Y,H=H,
                                                 kernel="Epanechnikov")$CV[,2])
      # warnings indicate Cholesky decomposition failed, in case when kNN is
      # smaller than the dimension of the model.
      # There, the CV criteria are are replaced by Inf
      CRIT.REG[Knn.grid<2] = Inf
      # here, the weight matrix had to be singular for kernels supported in
      # the interval [0,1]
      # plot(CRIT.REG~Knn.grid,type="b")
      if(rank.min != order(CRIT.REG)[1]) break
    }
    knn.reg = (Knn.grid)[order(CRIT.REG)[1]]
    #
    # knn for derivative
    #
    kNN.iterate = TRUE # indicator stopping the iteration in kNN CV procedure
    knn.max = 2        # initiate the kNN CV procedure
    ratio = 1
    while(kNN.iterate){
      knn.min = min(max(round(nlearn * Knn.learn.ratio[ratio]), knn.max),
                    nlearn-1)
      knn.max = min(max(round(nlearn * Knn.learn.ratio[ratio + 1]),
                        knn.min + kNN.grid.length  - 1),nlearn-1)
      ratio = ratio + 1
      if(ratio == length(Knn.learn.ratio)-1) kNN.iterate= FALSE # stop iteration
      if(knn.max == nlearn-1) kNN.iterate= FALSE                # stop iteration
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
      ##########################################################################
      # estimate regression operator with bandwidth maximizing the CV criterion
      ##########################################################################
      CRIT.REG = rep(NA,actual.kNN.grid.length)
      H = matrix(nrow=nlearn,ncol=actual.kNN.grid.length)
      for(knni in 1:actual.kNN.grid.length) H[,knni] = DIST.SORTED.ROW.BY.ROW[,
                                                             Knn.grid[knni] + 1]
      CRIT.REG = suppressWarnings(lqsm_cv_single(X,DIST,Y,H=H,
                                                 kernel="Epanechnikov")$CV[,2])
      # warnings indicate Cholesky decomposition failed, in case when kNN is
      # smaller than the dimension of the model.
      # There, the CV criteria are are replaced by Inf
      CRIT.REG[Knn.grid<3] = Inf
      # here, the weight matrix had to be singular for kernels supported in
      # the interval [0,1]
      # plot(CRIT.REG~Knn.grid,type="b")
      if(rank.min != order(CRIT.REG)[1]) break
    }
    knn.deriv = (Knn.grid)[order(CRIT.REG)[1]]
    ### CV choice of knn for derivative
    # results
    DIST.PRED.BY.LEARN = abs(outer(c(Xnew),c(X),"-"))
    DIST.SORTED.PRED.BY.LEARN = t(apply(DIST.PRED.BY.LEARN, 1, sort))
    res = llsm(X,Xnew,DIST.PRED.BY.LEARN,Y,
               h = DIST.SORTED.PRED.BY.LEARN[, knn.reg + 1], kernel=
                 "Epanechnikov")
    f[,k] = res[,1] # estimated responses
    g1[,k] = res[,2] # estimated derivatives based on local linear regression
    res.deriv = lqsm(X,Xnew,DIST.PRED.BY.LEARN,Y,
                     h = DIST.SORTED.PRED.BY.LEARN[, knn.deriv + 1], kernel=
                       "Epanechnikov")
    g[,k] = res.deriv[,2]
    ### visualisation of the fit
    # plot(Y~X)
    # points(res[,1]~Xnew,col=2,pch=16,cex=.5)
    # bw.CV0 <- regCVBwSelC(X, Y, deg=1, kernel=EpaK, interval=c(0,
    # diff(range(X))))
    # fit = locpol(Y~X,data = data.frame(X=X,Y = Y), deg=1, xeval=Xnew, kernel =
    # EpaK, bw=bw.CV0)
    # points(fit$lpFit[rank(Xnew),2]~Xnew, col=3,pch=16,cex=.5)
    # ## quadratic fitting
    # points(res.deriv[,1]~Xnew,col=4,pch=16,cex=.5)
    # bw.CV <- regCVBwSelC(X, Y, deg=2, kernel=EpaK, interval=c(0,
    # diff(range(X))))
    # fit = locpol(Y~X,data = data.frame(X=X,Y = Y), deg=2, xeval=Xnew, kernel =
    # EpaK, bw=bw.CV)
    # points(fit$lpFit[rank(Xnew),2]~Xnew, col=5, pch=16, cex=.5)
    # # ### derivative estimation
    # plot(res[,2]~Xnew)
    # points(fit$lpFit[rank(Xnew),3]~Xnew, col=3, pch=16, cex=.5)
    # points(res.deriv[,2]~Xnew,col=4,pch=16,cex=.5)
    # plot(Y-mean(Y)~ksi[,k])
    # lines(fit$lpFit[,2]~sort(Cnew[,k]))
  }
  return(list(reg.coefficients=f, deriv.coefficients=g, deriv.coefficients.lin =
                g1))
}


MuellerYao = function(C, Cnew, Y){
  # additive functional derivative estimation (Mueller and Yao ,2010)
  # library(locpol)
  ksi = C # t(t(C) - c(project(colMeans(X),phi)))
  nbasis = ncol(C)
  g = f = matrix(NA,nrow(Cnew),nbasis)
  bwd = rep(NA,nbasis)
  for(k in 1:nbasis){
    # regression estimation
    bw.CV0 <- regCVBwSelC(ksi[,k], Y-mean(Y), deg=1, kernel=EpaK, interval=c(0,
                                                          diff(range(ksi[,k]))))
    fit = locpol(Yc~ksi,data = data.frame(ksi=ksi[,k],Yc = Y-mean(Y)), deg=1,
                 xeval=Cnew[,k], kernel = EpaK, bw=bw.CV0)
    f[,k] = fit$lpFit[rank(Cnew[,k]),2]
    # derivative estimation
    bw.CV <- regCVBwSelC(ksi[,k], Y-mean(Y), deg=2, kernel=EpaK, interval=c(0,
                                                          diff(range(ksi[,k]))))
    fit = locpol(Yc~ksi,data = data.frame(ksi=ksi[,k],Yc = Y-mean(Y)), deg=2,
                 xeval=Cnew[,k], kernel = EpaK, bw=bw.CV)
    g[,k] = fit$lpFit[rank(Cnew[,k]),3]
    bwd[k] = bw.CV
  }
  g[is.na(g)] = 0
  f[is.na(f)] = 0
  return(list(coefficients=g,bwd=bwd,reg.coefficients=f))
}
