setwd("D:/Desktop/SF_LLR")
source("fllqr_funcs.R")

nsr = 0.40                               # noise-to-signal ratio
nlearn = 100                             # learning sample size
ntest = 500                              # testing sample size
Jmodel = 4                               # number of basis elements used for building functional predictors
################
dat = generate(nlearn, ntest, Jmodel, nsr)
LEARN = dat$LEARN
PRED = dat$PRED
Responses = dat$Responses
oracle.Responses = dat$oracle.Responses
oracle.DERIV = dat$oracle.DERIV
gridsize = ncol(LEARN)
BASIS = dat$BASIS
Regression = dat$fRegression
DERIV = dat$fDERIV
################
start.time <- Sys.time()
res1 = fllqr(Responses, LEARN, PRED)
end.time <- Sys.time()
res2 = fllr(Responses, LEARN, PRED)
## Plots ----
xrange = range(c(Regression[-(1:nlearn)]))
yrange = range(c(res1$Predicted.responses, res2$Predicted.responses))
par(mfrow = c(1, 1), mar = c(5, 5.35, 1, 1), cex.lab = 2, cex.axis = 1.5)
plot(
  Regression[-(1:nlearn)], res2$Predicted.responses, pch = 16, col = "grey", 
  xlim = xrange, ylim = yrange, 
  xlab = expression(m(X[i])), ylab = expression(hat(m)(X[i]))
)
points(Regression[-(1:nlearn)], res1$Predicted.responses, pch = 16, col = "red")
legend(
  "bottomright", c("Mean Regression", "Quantile Regression"), 
  pch = 16, col = c("grey", "red"), cex = 1.5
)

Grid = seq(0, 1, length = gridsize)
func_bounds <- range(cbind(res1$PREDICTED.DERIV[1:12, ], res2$PREDICTED.DERIV[1:12, ]))
par(mfrow = c(3,4), mar = c(2, 3, 1, 1), oma = c(0, 0, 0, 0), cex.axis = 1.5)
for(i in 1:12){
  plot(
    Grid, DERIV[i + nlearn, ], type = "l", lwd = 2, 
    ylim = func_bounds, bty = "l", ann = F
  )
  lines(Grid, res2$PREDICTED.DERIV[i, ], lwd = 3, lty = 2, col = "grey")
  lines(Grid, res1$PREDICTED.DERIV[i, ], lwd = 3, lty = 2, col = "red")
}

par(mfrow = c(1, 1))
############
difftime(end.time, start.time, units = "secs")



