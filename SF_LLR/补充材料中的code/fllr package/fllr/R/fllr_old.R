fllr_old = function(Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30, percent = 0.95,
                boot.seed = NULL, BASIS = NULL, method = "cv", derivative = TRUE, J.est = TRUE)
{
  # Check whether there are additional observations for predicting corresponding responses
  outsample = T
  if(identical(LEARN,PRED)) outsample = F
  ## if LEARN is the same matrix as PRED, outsample = F and PRED = LEARN
  ##
  # Estimate the second order moment operator
  gridsize = ncol(LEARN)
  if(is.null(BASIS)){
    if(outsample){
      sample.size = nrow(LEARN) + nrow(PRED)
      COV = crossprod(rbind(LEARN, PRED)) / (sample.size * (gridsize - 1))
    }else{
      sample.size = nrow(LEARN)
      COV = crossprod(LEARN) / (sample.size * (gridsize - 1))
    }
  }
  # Install the r package "parallelDist" to use "parDist" and compute
  # distances ||X_i - X_j||^2 between functional predictors
  require(parallelDist)
  require(fllr)
  DIST = as.matrix(parDist(LEARN, method = "euclidean")) / sqrt(gridsize - 1)
  DIST.SORTED.ROW.BY.ROW = t(apply(DIST, 1, sort))
  if(is.null(BASIS)){
    # Compute estimated first "dim.max" eigenfunctions
    res.eig = eigen(COV, sym = T)
    # Drop first eigenvalue and derive the approximating subspace expressing at least "percent"
    # of the remaining total variance
    cum.part.of.variance = cumsum(res.eig$values[-1]) / sum(res.eig$values[-1])
    dim.min = 1
    dim.max = sum(cum.part.of.variance < percent) + 2
    BASIS = t(res.eig$vectors[, 1:dim.max]) * sqrt(gridsize - 1)
  } else {
    # if J is estimated, dimensions are searched from one
    if(J.est) dim.min = 1 else dim.min = nrow(BASIS)
    dim.max = nrow(BASIS)
  }
  nlearn = nrow(LEARN)
  rank.min = kNN.grid.length
  Knn.learn.ratio = c(0.02, 0.25, 0.5, 0.75, 0.9)
  # Build an adaptative grid of kNN (number of knearest neighbours used for determining local bandwidths)
  # If the minimum of the criterion ("CRIT.REG") is reached for the bandwidth corresponding to
  # "Knn.grid[kNN.grid.length]", then the grid "Knn.grid" is updated
  for(ratio in 1:(length(Knn.learn.ratio) - 1)){
    knn.min = max(round(nlearn * Knn.learn.ratio[ratio]), dim.max + 4)
    knn.max = round(nlearn * Knn.learn.ratio[ratio + 1])
    knn.range = knn.max - knn.min
    if(knn.range <= kNN.grid.length){
      Knn.grid = knn.min:knn.max
    }else{
      step = ceiling(knn.range / kNN.grid.length)
      Knn.grid = seq(from = knn.min, to = knn.max, by = step)
    }
    kNN.grid.length = length(Knn.grid)
    rank.min = kNN.grid.length
    #########################################################################################
    # estimate regression operator with bandwidth maximizing cross validation criterion
    #########################################################################################
    CRIT.REG = matrix(Inf, kNN.grid.length, dim.max)
    dim.max = nrow(BASIS)
    B = BASIS[1:dim.max,]
    INNER.LEARN = tcrossprod(LEARN, B) / (gridsize - 1)
    if(method=="cv"){
      H = matrix(nrow=nrow(INNER.LEARN),ncol=length(Knn.grid))
      for(knni in 1:length(Knn.grid)) H[,knni] = DIST.SORTED.ROW.BY.ROW[, Knn.grid[knni] + 1]
      for(dim in dim.min:dim.max){
        CRIT.REG[,dim] = llsm_cv_single(matrix(INNER.LEARN[,1:dim],ncol=dim),DIST,Responses,H=H,kernel="Epanechnikov")$CV[,dim+1]
        dim.opt = dim
        if(dim >= 2){
          if(min(CRIT.REG[, dim - 1]) < min(CRIT.REG[, dim])){
            dim.opt = dim - 1
            break
          }
        }
      }
    }
    if(method=="aicc"){
      for(dim in (dim.min:dim.max)){
        # compute inner products between functional predictors and basis functions
        if(dim==1){
          B = as.matrix(t(BASIS[1, ]))
        }else{
          B = BASIS[1:dim, ]
        }
        INNER.LEARN = tcrossprod(LEARN, B) / (gridsize - 1)
        Estimated.responses = 0
        count = 0
        for(knn in Knn.grid){
          count = count + 1
          trHatMat = 0
          for(ii in 1:nlearn){
            PHI = cbind(1, t(t(INNER.LEARN) - INNER.LEARN[ii, ]))
            bwd = DIST.SORTED.ROW.BY.ROW[ii, knn + 1]
            U = DIST[ii, 1:nlearn] / bwd
            Kernel = (1 - U^2) * (U < 1)
            PHIK = PHI * Kernel
            TEMP = crossprod(PHIK, PHI)
            TEMP.INV = solve(TEMP)
            Temp = crossprod(PHIK, Responses)
            Estimated.responses[ii] = crossprod(TEMP.INV[1, ], Temp)
            # Compute the trace of the hat matrix for the AICc criterion
            trHatMat = trHatMat + crossprod(TEMP.INV[1, ], t(PHIK)[, ii])
          }
          # Compute criterion for selecting the optimal bandwidth
          Estimated.errors = Estimated.responses - Responses
          CRIT.REG[count, dim] = mean(Estimated.errors^2)
          CRIT.REG[count, dim] = log(CRIT.REG[count, dim]) + 1 + (2 * trHatMat + 2)/(nlearn - trHatMat - 2)
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
  mse = min(CRIT.REG[, dim.opt])
  # Determine optimal number of k-nearest neighbours and dimension J for the regression
  knn.reg = (Knn.grid)[order(CRIT.REG[, dim.opt])[1]]
  # Compute inner products between functional predictors and basis functions
  B = matrix(BASIS[1:dim.opt,],nrow=dim.opt)
  INNER.LEARN = matrix(INNER.LEARN[,1:dim.opt],ncol=dim.opt)
  if(method=="cv") Estimated.responses = llsm_leave(INNER.LEARN,INNER.LEARN,DIST,Responses,h = DIST.SORTED.ROW.BY.ROW[, knn.reg + 1], kernel="Epanechnikov")[,1]
  if(method=="aicc") Estimated.responses = llsm(INNER.LEARN,INNER.LEARN,DIST,Responses,h = DIST.SORTED.ROW.BY.ROW[, knn.reg + 1], kernel="Epanechnikov")[,1]
  ###########################################################################
  # estimate fuctional derivatives (with specific bandwidth choice)
  ###########################################################################
  if(derivative){
    # Step 1: compute estimated errors
    Estimated.errors = Responses - Estimated.responses
    # Set values for the wild bootstrap procedure
    gn = 0.5 * (1 + sqrt(5))
    gn.conj = 0.5 * (1 - sqrt(5))
    prob = (5 + sqrt(5)) / 10
    # Initialization
    Crit.deriv.boot = 0
    COEF.DERIV.BOOT = matrix(0, nlearn, dim.opt)
    PHIK = list()
    TEMP.INV = list()
    # Compute invariant quantities
    for(ii in 1:nlearn){
      PHI = cbind(1, t(t(INNER.LEARN) - INNER.LEARN[ii, ]))
      bwd = DIST.SORTED.ROW.BY.ROW[ii, knn.reg + 1]
      U = DIST[ii, 1:nlearn] / bwd
      Kernel = (1 - U^2) * (U < 1)
      Kernel[ii] = 0
      PHIK[[ii]] = PHI * Kernel
      TEMP = crossprod(PHIK[[ii]], PHI)
      TEMP.INV[[ii]] = solve(TEMP)
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
        COEF.DERIV.BOOT[ii, ] = COEF.DERIV.BOOT[ii, ] + TEMP.INV[[ii]][-1, ] %*% Temp
      }
    }
    # Compute the average of the boostrapped estimations for each functional predictor in the learning sample
    COEF.DERIV.BOOT = COEF.DERIV.BOOT / nboot
    # Compute the bootstrap pilot
    PILOT.BOOT.DERIV = COEF.DERIV.BOOT %*% B
    # Determine the optimal bandwidth for estimating the functional derivatives
    rank.min = kNN.grid.length
    Knn.learn.ratio = c(0.02, 0.25, 0.5, 0.75, 0.9)
    # Build an adaptative grid of kNN (number of knearest neighbours used for determining local bandwidths)
    # If the minimum of the criterion is reached for the last bandwidth
    # then the grid "Knn.grid" is updated
    for(ratio in 1:(length(Knn.learn.ratio) - 1)){
      knn.min = max(round(nlearn * Knn.learn.ratio[ratio]), dim.max + 4)
      knn.max = round(nlearn * Knn.learn.ratio[ratio + 1])
      knn.range = knn.max - knn.min
      if(knn.range <= kNN.grid.length){
        Knn.grid = knn.min:knn.max
      }else{
        step = ceiling(knn.range / kNN.grid.length)
        Knn.grid = seq(from = knn.min, to = knn.max, by = step)
      }
      kNN.grid.length = length(Knn.grid)
      rank.min = kNN.grid.length
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
  } else { # if the badnwidth for the derivative is just that from regression
    knn.deriv = knn.reg
  }
  # Save optimal bandwidths in terms of k-nearest neighbours into outputs
  knn.opt = list(reg = knn.reg, deriv = knn.deriv)
  COEF.DERIV = llsm(INNER.LEARN,INNER.LEARN,DIST,Responses,h=DIST.SORTED.ROW.BY.ROW[, knn.deriv + 1],kernel="Epa")[,-1]
  ESTIMATED.DERIV = COEF.DERIV %*% B
  # Compute predicted responses + predicted fuctional derivative if "outsample" is true
  # (i.e. one has additional observations for the functional predictor)
  if(outsample){
    DIST = as.matrix(parDist(rbind(LEARN, PRED), method = "euclidean")) / sqrt(gridsize - 1)
    DIST.SORTED.PRED.BY.LEARN = t(apply(DIST[-(1:nlearn), 1:nlearn], 1, sort))
    # Compute inner products between functional predictors and basis functions
    INNER.PRED = tcrossprod(PRED, B) / (gridsize - 1)
    nout = nrow(PRED)
    Predicted.responses = llsm(INNER.LEARN,INNER.PRED,DIST[-(1:nlearn), 1:nlearn],Responses,h=DIST.SORTED.PRED.BY.LEARN[,knn.reg + 1],kernel="Epa")[,1]
    COEF.PRED.DERIV = llsm(INNER.LEARN,INNER.PRED,DIST[-(1:nlearn), 1:nlearn],Responses,h=DIST.SORTED.PRED.BY.LEARN[,knn.deriv + 1],kernel="Epa")[,-1]
    PREDICTED.DERIV = COEF.PRED.DERIV %*% B
    return(list(Estimated.responses = Estimated.responses, ESTIMATED.DERIV = ESTIMATED.DERIV, Predicted.responses = Predicted.responses, PREDICTED.DERIV = PREDICTED.DERIV, mse = mse, knn.opt = knn.opt, dim.opt = dim.opt))
  }else{
    return(list(Estimated.responses = Estimated.responses, ESTIMATED.DERIV = ESTIMATED.DERIV, mse = mse, knn.opt = knn.opt, dim.opt = dim.opt))
  }
}
