library(fllr)
library(quantreg)
nsr = 0.4
gridsize = 101                           # grid size
Grid = seq(0, 1, length = gridsize)      # grid of measurements
nlearn = 100                             # learning sample size
ntest = 500                               # testing sample size
sample.size = nlearn + ntest             # whole sample size
Jmodel = 4                               # number of basis elements used for building functional predictors
BASIS = phif(Jmodel, Grid)               # first Jmodel Fourier basis elements
UNIF = matrix(runif(sample.size * Jmodel, min = -1, max = 1), nrow = sample.size, ncol = Jmodel)
PREDICTORS = UNIF %*% BASIS[1:Jmodel, ]

##########################################
# Simulate the regression and functional derivative
##########################################
TRANSFORMED.REG = exp(- UNIF^2)
TRANSFORMED.DERIV = - 2 * UNIF * TRANSFORMED.REG
Regression = apply(TRANSFORMED.REG, 1, sum)
DERIV = TRANSFORMED.DERIV %*% BASIS[1:Jmodel,]

###########################################
# Simulate scalar responses ("Responses")
###########################################
sigma = sqrt(nsr * var(Regression))
error = rnorm(sample.size, sd = sigma)
Responses.fullsample = Regression + error

###########################################
# Build learning sample
###########################################
LEARN = PREDICTORS[1:nlearn, ]
PRED = PREDICTORS[-(1:nlearn), ]
Responses = Responses.fullsample[1:nlearn]

nboot = 100
kNN.grid.length = 30
percent = 0.95
boot.seed = NULL
BASIS = NULL
method = "cv"
derivative = TRUE
J.est = TRUE


C = project(LEARN,BASIS)
Cnew = project(PRED,BASIS)
D = L2metric(PRED,LEARN)
res2 = llsm(C, Cnew, D, Responses, kNN = .05, kernel = "Gaussian")

h = apply(D, 1, function(x) quantile(x, 0.05))
res <- matrix(0, ntest, Jmodel + 1)
for(i in 1:ntest){
  temp <- t(t(C)-Cnew[i, ])
  kern_weight <- dnorm(D[i, ] / h[i])
  res[i, ] <- lm(Responses ~ temp, weights = kern_weight)$coef
}


D2 <- L2metric(LEARN,LEARN)
res2 <- llsm_leave(C, C, D2, Responses, kNN = .05, kernel = "Gaussian")
h = apply(D2, 1, function(x) quantile(x, 0.05))
res <- matrix(0, nlearn, Jmodel + 1)
for(i in 1:nlearn){
  temp <- t(t(C)-C[i, ])
  kern_weight <- dnorm(D2[i, ] / h[i])
  kern_weight[i] <- 0
  res[i, ] <- lm(Responses ~ temp, weights = kern_weight)$coef
}



res = fllr(Responses, LEARN, PRED)
res$knn.opt

plot(
  Regression[-(1:nlearn)], res$Predicted.responses, 
  xlab = expression(m(X[i])), ylab = expression(hat(m)(X[i])), 
  xlim = range(Regression[-(1:nlearn)]), ylim = range(res$Predicted.responses)
)


plot(DERIV[nlearn + 8, ], type = "l", lwd = 2)
lines(res$PREDICTED.DERIV[8, ], lwd = 2, lty = 2)





