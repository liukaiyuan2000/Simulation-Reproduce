################################################################################
# 
# Ferraty, F. and Nagy, S. (2021). Scalar-on-function local linear regression 
#                                  and beyond.
#
# Real data analysis: BERKELEY GROWTH DATASET (AOS HMY 2009)
#         Local linear estimation in the benchmark growth curves dataset 
#         and the estimation of the first functional derivatives of the 
#         regression function
#
################################################################################
library(fda)
library(fllr)
path2 = "graphics/"
# "growth$hgtm" = 31 by 39 numeric matrix giving the heights in centimeters of 
#                 39 boys at 31 ages
# "growth$hgtf" =  31 by 54 numeric matrix giving the heights in centimeters of 
#                  54 girls at 31 ages
# "growth$age" = numeric vector of length 31 giving the ages at which the 
#                heights were measured
Responses = growth$hgtm[31, ]
### put data into list format for computing growth velocity with "llderivgrowth"
GH = growth$hgtm[1:15, ]
nboys = ncol(GH)
exp.growthcurves = list()
exp.times = list()
for(nn in 1:nboys){
	exp.times[[nn]] = growth$age[1:15]
	exp.growthcurves[[nn]] = GH[, nn]
}
################################################################################
# smooth experimental growth profiles and their derivatives 
# (growth velocity profiles)
################################################################################
llderivgrowth = function(obs.growthcurves, obs.times, step = 1)
{
  ##############################################################################
  # "llderivgrowth" estimates growth curves and their 1st derivative by using 
  #                 the local linear regression 
  # "Responses" vector of length n containing the scalar responses Y1,...,Yn 
  # "obs.growthcurves" list containing k observed growth data: 
  #          obs.growthcurves[[1]]: X_1(t_{1,1}),...,X_1(t_{1,p_1})
  #          obs.growthcurves[[2]]: X_2(t_{2,1}),...,X_2(t_{2,p_2})
  #          .................................................
  #          obs.growthcurves[[k]]: X_k(t_{k,1}),...,X_k(t_{k,p_k})
  # "obs.times" list containing the grid of measurements (times) for each 
  #             observed growth data:
  #          obs.times[[1]]: t_{1,1},...,t_{1,p_1}
  #          obs.times[[2]]: t_{2,1},...,t_{2,p_2}
  #          .................................
  #          obs.times[[k]]: t_{k,1},...,t_{k,p_k}
  # "step" is the time unit allowing to build an equally-spaced grid points over 
  #          which the growth curves and their 1st derivatives are to be 
  #          estimated; "step" must match the real time unit used for the 
  #          varying temperature sequence (by default, "step" is set to 1)
  ##############################################################################
  # Returns a list containing
  # "$Eval" gives the equally-spaced grid points over which the 1st derivatives 
  #         are to be estimated
  # "$estimated.deriv.growth": list containing the estimated 1st derivatives of 
  #                            growth curves
  # "$estimated.growth": list containing the estimated growth curves
  ##############################################################################
  ## require "KernSmooth" and "sm" packages 
  require(KernSmooth)
  require(sm)
  K = length(obs.growthcurves)
  estimated.growthcurves = list()
  estimated.deriv.growthcurves = list()
  Eval = list()
  for(kk in 1:K){		
    MinMax = range(obs.times[[kk]])
    Eval[[kk]] = seq(MinMax[1], MinMax[2] + step, by = step)
    Range.grid = range(Eval[[kk]])
    hcv = h.select(obs.times[[kk]], obs.growthcurves[[kk]])
    estimated.growthcurves[[kk]] = locpoly(obs.times[[kk]], 
                                           obs.growthcurves[[kk]], drv = 0L, degree = 1, kernel = "normal", 
                                           bandwidth = hcv, 
                                           gridsize = length(Eval[[kk]]), bwdisc = 25, range.x = Range.grid, 
                                           binned = FALSE, truncate = TRUE)$y
    estimated.deriv.growthcurves[[kk]] = locpoly(obs.times[[kk]], 
                                                 obs.growthcurves[[kk]], drv = 1L, degree = 1, kernel = "normal", 
                                                 bandwidth = hcv, gridsize = length(Eval[[kk]]), bwdisc = 25, 
                                                 range.x = Range.grid, binned = FALSE, truncate = TRUE)$y
  }
  return(list(Eval = Eval, 
              estimated.deriv.growthcurves = estimated.deriv.growthcurves, 
              estimated.growthcurves = estimated.growthcurves))
}

# smoothing step with automatic bandwidth choice
res.llderivgrowth = llderivgrowth(exp.growthcurves, exp.times, step = 0.1)
################################################################################
# Functional local linear regression based on growth velocity profiles 
# as predictors
################################################################################
# use an ad hoc color palette
colfunc = colorRampPalette(rev(c("red", "orange", "pink", "green", "blue", 
                                 "brown")))
Colors = colfunc(nboys)
# display growth velocity profiles
pdf(file = paste(path2, "growth_velocity_profiles.pdf", sep = ""))
par(mfrow=c(1,1), mar = c(5,5,0.5,0.5), cex.lab = 2, cex.axis = 1.5)
xlim = range(exp.times)
ylim = range(res.llderivgrowth$estimated.deriv.growthcurves)
kk = 1
plot(res.llderivgrowth$Eval[[kk]], 
     res.llderivgrowth$estimated.deriv.growthcurves[[kk]],
	   xlab = "age (years)", ylab = "growth velocity", type = "l", lwd = 1, 
     xlim = xlim, ylim = ylim, col = Colors[kk])
for(kk in 2:nboys)
	lines(res.llderivgrowth$Eval[[kk]], 
	      res.llderivgrowth$estimated.deriv.growthcurves[[kk]], lwd = 1, 
	      col = Colors[kk])
dev.off()
# put growth velocity profiles into a matrix
GVP = matrix(0, nboys, length(res.llderivgrowth$Eval[[1]]))
for(nn in 1:nboys) GVP[nn, ] = 
  res.llderivgrowth$estimated.deriv.growthcurves[[nn]]
res = fllr(Responses, GVP, GVP, nboot = 100, kNN.grid.length = 30, 
           percent = 0.95)
# use an ad hoc color palette
colfunc = colorRampPalette(rev(c("red", "orange", "pink", "green", "blue", 
                                 "brown")))
Colors = colfunc(nboys)
# display predicted Frechet derivatives
pdf(file = paste(path2, "growth_estimated_fderiv.pdf", sep = ""))
par(mfrow = c(1,1), mar = c(5,5,0.5,0.5), cex.axis = 1.5, cex.lab = 2)
ylim = range(res$ESTIMATED.DERIV)
plot(res.llderivgrowth$Eval[[1]], res$ESTIMATED.DERIV[1, ], xlab= "age (years)",
     ylab = "functional derivative", lwd = 0.5, col = Colors[1], type = "l", 
     ylim = ylim)
for(ii in 2:nboys) lines(res.llderivgrowth$Eval[[1]], res$ESTIMATED.DERIV[ii, ],
                         col = Colors[ii])
dev.off()
# Display estimated heigths against observed ones
pdf(file = paste(path2, "growth_obs_vs_estimations.pdf", sep = ""))
par(mfrow = c(1,1), mar = c(5,5,0.5,0.5), cex.axis = 1.5, cex.lab = 2)
plot(Responses, res$Estimated.responses, xlab = "height (cm) at 18",
     ylab = "Estimated height (cm) at 18")
abline(0,1)
dev.off()
################################################################################
# Single functional index regression model based on 1-10 (whole) growth profiles
# as predictors
################################################################################
# "growth$hgtm" = 31 by 39 numeric matrix giving the heights in centimeters of 
#                 39 boys at 31 ages
# "growth$hgtf" =  31 by 54 numeric matrix giving the heights in centimeters of 
#                  54 girls at 31 ages
# "growth$age" = numeric vector of length 31 giving the ages at which the 
#                heights were measured
Responses = growth$hgtm[31, ]
### put data into list format for computing growth velocity with "llderivgrowth"
GH = growth$hgtm[1:15, ]
nboys = ncol(GH)
exp.growthcurves = list()
exp.times = list()
for(nn in 1:nboys){
	exp.times[[nn]] = growth$age[1:15]
	exp.growthcurves[[nn]] = GH[, nn]
}
# smooth experimental growth profiles and their derivatives 
# (growth velocity profiles)
# smoothing step with automatic bandwidth choice
res.llderivgrowth = llderivgrowth(exp.growthcurves, exp.times, step = 0.1)
# put growth velocity profiles into a matrix
GVP = matrix(0, nboys, length(res.llderivgrowth$Eval[[1]]))
for(nn in 1:nboys) GVP[nn, ] = 
  res.llderivgrowth$estimated.deriv.growthcurves[[nn]]
# local linear estimation of responses and  functional derivatives
res = fllr(Responses, GVP, GVP, nboot = 100, kNN.grid.length = 30, 
           percent = 0.95)
# single functional index model
FuncIndexEst = colSums(res$ESTIMATED.DERIV) / nboys
gridsize = length(res.llderivgrowth$Eval[[1]])
FuncIndexEst = FuncIndexEst / sqrt((sum(FuncIndexEst^2) - 
                              0.5 * sum(FuncIndexEst[c(1, gridsize)]^2)) / 0.1)
Proj.coef = (GVP %*% FuncIndexEst - 
           0.5 * GVP[, c(1, gridsize)] %*% FuncIndexEst[c(1, gridsize)]) / 0.1
# require "KernSmooth" package
require(KernSmooth)
# regress "Responses" on "Proj.coef"
hcv = dpill(Proj.coef, Responses)
reslocpoly = locpoly(Proj.coef, Responses, drv = 0L, degree = 1, 
                     kernel = "normal", bandwidth = hcv,
                     gridsize = 1000, bwdisc = 25, range.x = range(Proj.coef), 
                     binned = FALSE, truncate = TRUE)
Estimated.responses = 0
for(jj in 1:nboys){
	ind.nearest.values.in.x = which.min(abs(Proj.coef[jj] - reslocpoly$x))
	Estimated.responses[jj] = reslocpoly$y[ind.nearest.values.in.x]
}
# Display estimated heights against observed ones
pdf(file = paste(path2, "growth_sfim_obs_vs_estimations.pdf", sep = ""))
par(mfrow = c(1,1), mar = c(5,5,0.5,0.5), cex.lab = 2, cex.axis = 1.5)
df = data.frame(x = Responses, y = Estimated.responses)
with(df, symbols(x = x, y=y, circles = rep(3, nboys), inches = 1/20, 
                 bg = "black", xlab = "height (cm) at 18", 
                 ylab = "Estimated height (cm) at 18"))
abline(0,1)
dev.off()
# Display estimated link function
pdf(file = paste(path2, "growth_sfim_linkfunction.pdf", sep = ""))
par(mfrow = c(1,1), mar = c(5,3,0.5,0.5), cex.lab = 2, cex.axis = 1.5)
Bounds.y = range(c(reslocpoly$y, Responses))
plot(Proj.coef, Responses, xlab = "", ylab = "Estimated link function", 
     ylim = c(min(Proj.coef), 200))
lines(reslocpoly, lwd = 5)
dev.off()
# Display functional index in "blue"
pdf(file = paste(path2, "growth_sfim_functionalindex.pdf", sep = ""))
par(mfrow = c(1,1), mar = c(5,3,0.5,0.5), cex.lab = 2, cex.axis = 1.5)
Bounds.y = range(FuncIndexEst)
plot(res.llderivgrowth$Eval[[1]], FuncIndexEst, xlab = "age (years)", ylab = "",
     ylim = Bounds.y, type = "l", lwd = 5, col = "blue")
dev.off()
# Display functional index in "black"
pdf(file = paste(path2, "growth_sfim_functionalindex_black.pdf", sep = ""))
par(mfrow = c(1,1), mar = c(5,3,0.5,0.5), cex.lab = 2, cex.axis = 1.5)
Bounds.y = range(FuncIndexEst)
plot(res.llderivgrowth$Eval[[1]], FuncIndexEst, xlab = "age (years)", ylab = "", 
     ylim = Bounds.y, type = "l", lwd = 5, col = "black")
dev.off()
# estimation based on 6-10 growth velocity
GVP.AFTER.6 = GVP[, -(1:50)]
# local linear estimation of responses and  functional derivatives
res.after.6 = fllr(Responses, GVP.AFTER.6, GVP.AFTER.6, nboot = 100, 
                   kNN.grid.length = 30, percent = 0.95)
# single functional index model
FuncIndexEst.after.6 = colSums(res.after.6$ESTIMATED.DERIV) / nboys
gridsize.after.6 = ncol(res.after.6$ESTIMATED.DERIV)
FuncIndexEst.after.6 = FuncIndexEst.after.6 / sqrt((sum(FuncIndexEst.after.6^2)-
             0.5 * sum(FuncIndexEst.after.6[c(1, gridsize.after.6)]^2)) / 0.1)
Proj.coef.after.6 = (GVP.AFTER.6 %*% FuncIndexEst.after.6 - 
                       0.5 * GVP.AFTER.6[, c(1, gridsize.after.6)] %*% 
                       FuncIndexEst.after.6[c(1, gridsize.after.6)]) / 0.1
hcv = dpill(Proj.coef.after.6, Responses)
reslocpoly = locpoly(Proj.coef.after.6, Responses, drv = 0L, degree = 1, 
                     kernel = "normal", bandwidth = hcv,
                     gridsize = 1000, bwdisc = 25, 
                     range.x = range(Proj.coef.after.6), binned = FALSE, 
                     truncate = TRUE)
Estimated.responses.after.6 = 0
for(jj in 1:nboys){
	ind.nearest.values.in.x = which.min(abs(Proj.coef.after.6[jj] - reslocpoly$x))
	Estimated.responses.after.6[jj] = reslocpoly$y[ind.nearest.values.in.x]
}
cor(Estimated.responses.after.6, Responses)
# compare estimation based on 1-10 growth velocity or 6-10 growth velocity
pdf(file = paste(path2, 
  "growth_sfim_obs_vs_estimations_compare_1-10_6-10_white_black.pdf", sep = ""))
Ybounds = range(Estimated.responses, Estimated.responses.after.6)
par(mfrow = c(1,1), mar = c(5,5,1,0.5), cex.lab = 2, cex.axis = 1.5)
df = data.frame(x = Responses, y = Estimated.responses.after.6)
with(df, symbols(x = x, y=y, circles = rep(3, nboys), inches = 1/18, 
                 bg = "white", fg = "black", xlab = "height (cm) at 18", 
                 ylab = "Estimated height (cm) at 18", ylim = Ybounds))
abline(0,1)
points(Responses, Estimated.responses, pch = 19, col = "black")
dev.off()
# compute l2 norm of the predicted Frechet derivative
gridsize = ncol(res$ESTIMATED.DERIV)
L2norm.est = sqrt((apply(res$ESTIMATED.DERIV^2, 1, sum) - 
             0.5 * apply(res$ESTIMATED.DERIV[, c(1, gridsize)]^2, 1, sum))/ 0.1)