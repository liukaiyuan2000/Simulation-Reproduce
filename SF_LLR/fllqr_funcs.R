rm(list = ls())
library(quantreg)
library(fllr)
library(ggplot2)
Rcpp::sourceCpp("D:/Desktop/学习/SF_LLR/fllqr_funcs.cpp")

fllqr <- function(
    Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30,
    percent = 0.95, tau = 0.5
  ){
  gridsize = ncol(LEARN)
  sample.size = nrow(LEARN) + nrow(PRED)
  ## Cov/100
  COV = crossprod(rbind(LEARN, PRED))/(sample.size * (gridsize - 1))

  require(parallelDist)
  require(fllr)
  DIST = as.matrix(parDist(LEARN, method = "euclidean"))/sqrt(gridsize - 1)
  DIST.SORTED.ROW.BY.ROW = t(apply(DIST, 1, sort))

  res.eig = eigen(COV, sym = T)
  cum.part.of.variance = cumsum(res.eig$values[-1])/sum(res.eig$values[-1])
  dim.min = 1
  dim.max = sum(cum.part.of.variance < percent) + 2
  BASIS = t(res.eig$vectors[, 1:dim.max]) * sqrt(gridsize - 1)

  nlearn = nrow(LEARN)
  rank.min = kNN.grid.length
  Knn.learn.ratio = c(0.02, 0.25, 0.5, 0.75, 0.9)
  for (ratio in 1:(length(Knn.learn.ratio) - 1)) {
    knn.min = max(round(nlearn * Knn.learn.ratio[ratio]), dim.max + 4)
    knn.max = round(nlearn * Knn.learn.ratio[ratio + 1])
    knn.range = knn.max - knn.min
    if (knn.range <= kNN.grid.length) {
      Knn.grid = knn.min:knn.max
    } else {
      step = ceiling(knn.range/kNN.grid.length)
      Knn.grid = seq(from = knn.min, to = knn.max, by = step)
    }
    kNN.grid.length = length(Knn.grid)
    rank.min = kNN.grid.length
    CRIT.REG = matrix(Inf, kNN.grid.length, dim.max)
    dim.max = nrow(BASIS)
    B = BASIS[1:dim.max, ]
    INNER.LEARN = tcrossprod(LEARN, B)/(gridsize - 1)
    H = matrix(nrow = nrow(INNER.LEARN), ncol = length(Knn.grid))
    for (knni in 1:length(Knn.grid)) {
      H[, knni] = DIST.SORTED.ROW.BY.ROW[, Knn.grid[knni] + 1]
    }
    for (dim in dim.min:dim.max) {
      CRIT.REG[, dim] = llqr_cv_cpp(
        cbind(1, matrix(INNER.LEARN[, 1:dim], ncol = dim)), DIST, Responses, H = H,
        kernI = 3, tau)[, dim + 1]
      dim.opt = dim
      if (dim >= 2) {
        if (min(CRIT.REG[, dim - 1]) < min(CRIT.REG[, dim])) {
          dim.opt = dim - 1
          break
        }
      }
    }

    if (rank.min != order(CRIT.REG[, dim.opt])[1])
      break
  }
  mse = min(CRIT.REG[, dim.opt])
  knn.reg = (Knn.grid)[order(CRIT.REG[, dim.opt])[1]]
  B = matrix(BASIS[1:dim.opt, ], nrow = dim.opt)
  INNER.LEARN = matrix(INNER.LEARN[, 1:dim.opt], ncol = dim.opt)
  Estimated.responses = llqr_leave_cpp(
    cbind(1, INNER.LEARN), cbind(1, INNER.LEARN), DIST, Responses,
    h = DIST.SORTED.ROW.BY.ROW[, knn.reg + 1], kernI = 3, tau
  )[, 1]

  Estimated.errors = Responses - Estimated.responses
  gn = 0.5 * (1 + sqrt(5))
  gn.conj = 0.5 * (1 - sqrt(5))
  prob = (5 + sqrt(5))/10
  Crit.deriv.boot = 0
  COEF.DERIV.BOOT = matrix(0, nlearn, dim.opt)
  PHIK = list()
  TEMP.INV = list()
  for (ii in 1:nlearn) {
    PHI = cbind(1, t(t(INNER.LEARN) - INNER.LEARN[ii,
    ]))
    bwd = DIST.SORTED.ROW.BY.ROW[ii, knn.reg + 1]
    U = DIST[ii, 1:nlearn]/bwd
    Kernel = (1 - U^2) * (U < 1)
    Kernel[ii] = 0
    PHIK[[ii]] = PHI * Kernel
    TEMP = crossprod(PHIK[[ii]], PHI)
    TEMP.INV[[ii]] = solve(TEMP)
  }
  for (bb in 1:nboot) {
    U = runif(nlearn)
    V = (U >= prob) * gn + (U < prob) * gn.conj
    Boostrapped.errors = Estimated.errors * V
    Boostrapped.responses = Estimated.responses + Boostrapped.errors
    for (ii in 1:nlearn) {
      Temp = crossprod(PHIK[[ii]], Boostrapped.responses)
      COEF.DERIV.BOOT[ii, ] = COEF.DERIV.BOOT[ii, ] +
        TEMP.INV[[ii]][-1, ] %*% Temp
    }
  }
  COEF.DERIV.BOOT = COEF.DERIV.BOOT/nboot
  PILOT.BOOT.DERIV = COEF.DERIV.BOOT %*% B
  rank.min = kNN.grid.length
  Knn.learn.ratio = c(0.02, 0.25, 0.5, 0.75, 0.9)
  for (ratio in 1:(length(Knn.learn.ratio) - 1)) {
    knn.min = max(round(nlearn * Knn.learn.ratio[ratio]), dim.max + 4)
    knn.max = round(nlearn * Knn.learn.ratio[ratio + 1])
    knn.range = knn.max - knn.min
    if (knn.range <= kNN.grid.length) {
      Knn.grid = knn.min:knn.max
    } else {
      step = ceiling(knn.range/kNN.grid.length)
      Knn.grid = seq(from = knn.min, to = knn.max, by = step)
    }
    kNN.grid.length = length(Knn.grid)
    rank.min = kNN.grid.length
    Crit.deriv.boot = rep(NA, kNN.grid.length)
    count = 0
    for (knn in Knn.grid) {
      count = count + 1
      COEF.DERIV = llqr_leave_cpp(
        cbind(1, INNER.LEARN), cbind(1, INNER.LEARN), DIST, Responses,
        h = DIST.SORTED.ROW.BY.ROW[, knn + 1], kernI = 3, tau
      )[, -1]
      ESTIMATED.DERIV = COEF.DERIV %*% B
      Crit.deriv.boot[count] = mean((PILOT.BOOT.DERIV - ESTIMATED.DERIV)^2)
    }
    if (rank.min != order(Crit.deriv.boot)[1])
      break
  }
  crit.deriv.boot = min(Crit.deriv.boot)
  knn.deriv = (Knn.grid)[which.min(Crit.deriv.boot)]
  knn.opt = list(reg = knn.reg, deriv = knn.deriv)
  COEF.DERIV = llqr_cpp(
    cbind(1, INNER.LEARN), cbind(1, INNER.LEARN), DIST, Responses,
    h = DIST.SORTED.ROW.BY.ROW[, knn.deriv + 1], kernI = 3, tau
  )[, -1]
  ESTIMATED.DERIV = COEF.DERIV %*% B

  DIST = as.matrix(parDist(rbind(LEARN, PRED), method = "euclidean"))/sqrt(gridsize - 1)
  DIST.SORTED.PRED.BY.LEARN = t(apply(DIST[-(1:nlearn), 1:nlearn], 1, sort))
  INNER.PRED = tcrossprod(PRED, B)/(gridsize - 1)
  nout = nrow(PRED)
  Predicted.responses = llqr_cpp(
    cbind(1, INNER.LEARN), cbind(1, INNER.PRED), DIST[-(1:nlearn), 1:nlearn], Responses,
    h = DIST.SORTED.PRED.BY.LEARN[, knn.reg + 1], kernI = 3, tau
  )[, 1]
  COEF.PRED.DERIV = llqr_cpp(
    cbind(1, INNER.LEARN), cbind(1, INNER.PRED), DIST[-(1:nlearn), 1:nlearn], Responses,
    h = DIST.SORTED.PRED.BY.LEARN[, knn.reg + 1], kernI = 3, tau
  )[, -1]
  PREDICTED.DERIV = COEF.PRED.DERIV %*% B
  return(
    list(
      Estimated.responses = Estimated.responses,
      ESTIMATED.DERIV = ESTIMATED.DERIV, Predicted.responses = Predicted.responses,
      PREDICTED.DERIV = PREDICTED.DERIV, mse = mse, knn.opt = knn.opt,
      dim.opt = dim.opt
    )
  )
}


#' Data simulation
#'
#' Simulates a random sample from models (M1)-(M3) from Ferraty and Nagy (2020).
#'
#' @param nlearn    number of learning functions
#' @param ntest     number of testing functions
#' @param Jmodel    effective dimension of the regressors, number of basis
#' functions on which the responses truly depend
#' @param nsr       noise-to-signal ration, a number between 0 and 1
#' @param aa        coefficient a in model \code{M3}, determining the linear
#' combination of the linear and the non-linear part of the regression
#' @param rho       parameter that drives the noise in the confounding part of
#' the regressors in models \code{M2} and \code{M3}, i.e. the part of the
#' regressors that corresponds to basis functions that do not contribute to the
#' modelling of the response
#' @param model     indicator of the model used, accepted are \code{"M1"},
#' \code{"M2"} or \code{"M3"}
#'
#' @return A list with the following components:
#' \itemize{
#' \item \code{LEARN} Matrix of dimension nlearn-times-101, one learning
#' function per row.
#' \item \code{PRED}  Matrix of new regressors of dimension npred-times-101, one
#' function per row.
#' \item \code{Responses} Vector of responses that correspond to the rows of
#' \code{LEARN}.
#' \item \code{oracle.Responses} True (unknown) responses generated by the model
#' that correspond to the rows of \code{PRED}.
#' \item \code{oracle.DERIV} True (unknown) functional derivatives generated
#' by the model that correspond to the rows of \code{PRED}.
#' \item \code{fDERIV = DERIV} The complete matrix of functional derivatives of
#' dimension (nlearn+npred)-times-101, first nlearn functions correspond to the
#' rows of \code{LEARN}, the remaining npred function to the rows of
#' \code{PRED}.
#' \item \code{fRegression} The complete vector of noiseless conditional
#' expectations in the model of length (nlearn+npred), first nlearn values
#' correspond to the rows of \code{LEARN}, the remaining npred values to the
#' rows of \code{PRED}.
#' \item \code{DERIV} True functional derivatives generated by the model.
#' Matrix of dimension nlearn-times-101, one row per a row of \code{LEARN}.
#' \item \code{BASIS} Matrix of the basis vectors into which the functional
#' regressors were decomposed.
#' }
#'
#' @author Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
#'
#'	@references Ferraty, F., and Nagy, S. (2020).
#'	Scalar-on-function local linear regression and beyond.
#'	\emph{Under review}.
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

generate = function(nlearn=100, ntest=500, Jmodel=4, nsr=.05, aa=1, rho = NULL,
                    model="M1"){
  # simulates a random sample in the study in the paper,
  # using the parameters provided
  # nlearn - learning sample size
  # ntest  - testing sample size
  # Jmodel - effective dimension of the dataset, Response depends only on Jmodel
  #          pricipal components of the data
  # nsr    - noise to signal ratio
  # aa     - parameter a tuning for non-linearity of the regression problem
  # model  - "M1" - J basis functions, only non-linear response
  #        - "M2" - 2*J basis functions, only non-linear response
  #        - "M3" - 2*J basis functions, linear and non-linear response

  gridsize = 101                           # grid size
  Grid = seq(0, 1, length = gridsize)      # grid of measurements
  sample.size = nlearn + ntest             # whole sample size
  if(model=="M1") BASIS = phif(Jmodel, Grid)
  if(model=="M2"|model=="M3") BASIS = phif(2 * Jmodel, Grid)
  # if model M2 and M3 compute 2*J Fourier basis elements
  bound = function(ratio){sqrt(ratio/(1 - ratio))}
  ################################
  # Simulate functional predictors
  ################################
  if(is.null(rho)) rho = nsr
  UNIF1 = matrix(runif(sample.size * Jmodel, min = -1, max = 1), nrow =
                   sample.size, ncol = Jmodel)
  b = bound(rho)
  UNIF2 = matrix(runif(sample.size * Jmodel, min = -b, max = b), nrow =
                   sample.size, ncol = Jmodel)
  if(model=="M1") PREDICTORS = UNIF1 %*% BASIS[1:Jmodel, ]
  if(model=="M2"|model=="M3") PREDICTORS = UNIF1 %*% BASIS[1:Jmodel, ] +
    UNIF2 %*% BASIS[-(1:Jmodel), ]
  ##############################################################################
  # Simulate linear part for regression and functional derivative
  ##############################################################################
  Beta = apply(BASIS[1:Jmodel, ], 2, sum)
  linear.part = (PREDICTORS %*% Beta - 0.5 * PREDICTORS[, c(1, gridsize)] %*%
                   Beta[c(1, gridsize)]) / (gridsize - 1)
  DERIV.LINEAR = matrix(rep(Beta, sample.size), nrow = sample.size, ncol =
                          gridsize, byrow = T)
  ##############################################################################
  # Simulate nonparametric part for regression and functional derivative
  ##############################################################################
  TRANSFORMED.REG = exp(- UNIF1^2)
  TRANSFORMED.DERIV = - 2 * UNIF1 * TRANSFORMED.REG
  np.part = apply(TRANSFORMED.REG, 1, sum)
  DERIV.NP = TRANSFORMED.DERIV %*% BASIS[1:Jmodel,]
  ##############################################################################
  # Compute the whole regression operator and functional derivaitve
  ##############################################################################
  if(model=="M1"|model=="M2"){
    Regression = np.part
    DERIV = DERIV.NP
  }
  if(model=="M3"){
    Regression = (1 - aa) * linear.part + aa * np.part
    DERIV = (1 - aa) * DERIV.LINEAR + aa * DERIV.NP
  }
  ##############################################################################
  # Simulate scalar responses (Responses) with given signal-to-signal ratio nsr
  ##############################################################################
  sigma = sqrt(nsr * var(Regression))
  error = rnorm(sample.size, sd = sigma)
  Responses.fullsample = Regression + error
  #####################################
  # Build learning sample
  #####################################
  LEARN = PREDICTORS[1:nlearn, ]
  PRED = PREDICTORS[-(1:nlearn), ]
  Responses = Responses.fullsample[1:nlearn]
  return(list(LEARN=LEARN, PRED = PRED, Responses = Responses,
              oracle.Responses = Regression[-(1:nlearn)],
              oracle.DERIV = DERIV[-(1:nlearn), ],
              fDERIV = DERIV,
              fRegression = Regression,
              DERIV = DERIV[(1:nlearn), ],
              BASIS = BASIS))
}

