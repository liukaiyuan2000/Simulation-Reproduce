set.seed(10)
res = list()
fRegression = list()
fDERIV = list()
fnsr = rev(range(Nsr.grid))
fnlearn = c(100,500)
for(i in 1:2){
  nsr = fnsr[i]
  nlearn = fnlearn[i]
  sample.size = nlearn + ntest
  ################################
  # Simulate all
  ################################
  dat = generate(nlearn, ntest, Jmodel, nsr, model="M1")
  LEARN = dat$LEARN
  PRED = dat$PRED
  Responses = dat$Responses
  oracle.Responses = dat$oracle.Responses
  oracle.DERIV = dat$oracle.DERIV
  gridsize = ncol(LEARN)
  BASIS = dat$BASIS
  ######################################
  # Variances for denominators of ORMSEP
  ######################################
  denom.reg = mean((oracle.Responses-mean(oracle.Responses))^2)
  denom.deriv = mean((t(oracle.DERIV) - colMeans(oracle.DERIV))^2)
  ##############################################
  # Functional local linear estimator
  ##############################################
  res[[i]] = fllr(Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30, 
                  percent = 0.95, BASIS = BASIS, boot.seed=1, J.est = FALSE)
  resB = fllr(Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30, 
              percent = 0.95, BASIS = BASIS, boot.seed = 1, derivative=FALSE, 
              J.est = FALSE)
  res[[i]]$PREDICTED.DERIV.HREG = resB$PREDICTED.DERIV
  #
  fRegression[[i]] = dat$fRegression
  fDERIV[[i]] = dat$fDERIV
}

xbounds = range(c(fRegression[[1]][-(1:(fnlearn[1]))],
                  fRegression[[2]][-(1:(fnlearn[2]))]))
ybounds = range(c(res[[1]]$Predicted.responses,res[[2]]$Predicted.responses))

par(mfrow = c(1, 1), mar = c(5, 5.35, 1, 1), cex.lab = 2, cex.axis = 1.5)
plot(fRegression[[1]][-(1:(fnlearn[1]))], res[[1]]$Predicted.responses, 
     xlab = expression(m(X[i])), ylab = expression(hat(m)(X[i])), 
     xlim = xbounds, ylim = ybounds)

par(mfrow = c(1, 1), mar = c(5, 5.35, 1, 1), cex.lab = 2, cex.axis = 1.5)
plot(fRegression[[2]][-(1:(fnlearn[2]))], res[[2]]$Predicted.responses, 
     xlab = expression(m(X[i])), ylab = expression(hat(m)(X[i])), 
     xlim = xbounds, ylim = ybounds)



Grid = seq(0, 1, length = gridsize)
fignames = c(paste0("deriv_predictions_hboot_vs_hreg_n",fnlearn[1],"_nsr04"),
             paste0("deriv_predictions_hboot_vs_hreg_n",fnlearn[2],"_nsr005"))
for(i in 1:2){
  pdf(file = paste0("graphics/",fignames[i],".pdf"))
  Sample.index = 1:12
  Bounds = range(cbind(res[[i]]$PREDICTED.DERIV[Sample.index, ], 
                       res[[i]]$PREDICTED.DERIV.HREG[Sample.index, ], 
                       fDERIV[[i]][Sample.index + fnlearn[i], ]))
  par(mfrow = c(3,4), mar = c(2, 3, 1, 1), oma = c(0, 0, 0, 0), cex.axis = 1.5)
  for(ii in Sample.index){
    plot(Grid, fDERIV[[i]][ii + fnlearn[i], ], type = "l", lwd = 2, 
         ylim = Bounds, bty = "l", ann = F)
    # lines(Grid, res[[i]]$PREDICTED.DERIV.HREG[ii, ], lwd = 3, lty = 3, 
    # col = "gray")
    lines(Grid, res[[i]]$PREDICTED.DERIV[ii, ], lwd = 3, lty = 2, col = "red")
  }
  dev.off()
}




