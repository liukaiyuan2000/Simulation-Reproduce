################################################################################
# 
# Ferraty, F. and Nagy, S. (2021). Scalar-on-function local linear regression 
#                                  and beyond.
#
#  Model (M1): Figures and a table with bandwidths
#  J=4 and the Fourier basis for building the regression model
#  Fourier basis is used to predict the functions
#  data-driven choice of 
#		* two bandwidths: - one for the regression operator and 
#                     - one for its functional derivatives
#
################################################################################

set.seed(1)
Jmodel = 4
nlearn = n.max                           # testing sample size
sample.size = nlearn + ntest             # whole sample size
nbsamples = NBsamples                    # number of independent runs
#########################################
#
# START LOOP
#
#########################################
Nsr.grid = c(0.05, 0.2, 0.4)
n.grid = seq(n.min, n.max, by=50) # sample sizes
TAB.RES.MED = array(NA, c(2, 8, length(n.grid), length(Nsr.grid)))
TAB.RES.MEAN = array(NA, c(2, 8, length(n.grid), length(Nsr.grid)))
##
count.n = 0
for(nlearn in n.grid){
	count.n = count.n + 1
	sample.size = nlearn + ntest
	count.nsr = 0
	for(nsr in Nsr.grid){
		count.nsr = count.nsr + 1
		cat("n = ", nlearn, " and nsr = ", nsr, "\n", sep ="")
		RESULTS = foreach(kk = 1:nbsamples, .combine='rbind') %dopar% {
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
			### cv
			res = fllr(Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30,
			           percent = 0.95, BASIS = BASIS, J.est = FALSE)
			msep.reg.oracle.cv = mean((oracle.Responses - res$Predicted.responses)^2)/
			  denom.reg
			msep.deriv.oracle.cv = mean((oracle.DERIV - res$PREDICTED.DERIV)^2)/
			  denom.deriv
			### aicc
			res = fllr(Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30,
			           percent = 0.95, BASIS = BASIS, method="aicc", J.est = FALSE)
			msep.reg.oracle.aicc = 
			  mean((oracle.Responses - res$Predicted.responses)^2)/denom.reg
			msep.deriv.oracle.aicc = mean((oracle.DERIV - res$PREDICTED.DERIV)^2)/
			  denom.deriv
			### without estimation of the bandwidth for derivative
			### cv
			res = fllr(Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30, 
			           percent = 0.95, BASIS = BASIS, derivative=FALSE, J.est = FALSE)
			msep.reg.oracle.cv.nd = 
			  mean((oracle.Responses - res$Predicted.responses)^2)/denom.reg
			msep.deriv.oracle.cv.nd = mean((oracle.DERIV - res$PREDICTED.DERIV)^2)/
			  denom.deriv
			### aicc
			res = fllr(Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30, 
			           percent = 0.95, BASIS = BASIS, method="aicc", derivative=FALSE,
			           J.est = FALSE)
			msep.reg.oracle.aicc.nd = 
			  mean((oracle.Responses - res$Predicted.responses)^2)/denom.reg
			msep.deriv.oracle.aicc.nd = mean((oracle.DERIV - res$PREDICTED.DERIV)^2)/
			  denom.deriv
			###
			return(c(msep.reg.oracle.cv, msep.deriv.oracle.cv,
			         msep.reg.oracle.aicc, msep.deriv.oracle.aicc,
			         msep.reg.oracle.cv.nd, msep.deriv.oracle.cv.nd,
			         msep.reg.oracle.aicc.nd, msep.deriv.oracle.aicc.nd))
		}
	TAB.RES.MED[, , count.n, count.nsr] = t(cbind(apply(RESULTS, 2, median), 
	                                              apply(RESULTS, 2, mad)))
	TAB.RES.MEAN[, , count.n, count.nsr] = t(cbind(apply(RESULTS, 2, mean), 
	                                               apply(RESULTS, 2, sd)))
}
}
#########################################
#
# END LOOP
#
#########################################
dimnames(TAB.RES.MEAN)[[1]] = c("mean", "sd")
dimnames(TAB.RES.MED)[[1]] = c("median", "mad")
dimnames(TAB.RES.MEAN)[[2]] = dimnames(TAB.RES.MED)[[2]] = c(
  "msep reg oracle cv", "msep deriv oracle cv", "msep reg oracle aicc", 
  "msep deriv oracle aicc", "msep reg oracle cv non-der", 
  "msep deriv oracle cv non-der", "msep reg oracle aicc non-der", 
  "msep deriv oracle aicc non-der")
dimnames(TAB.RES.MEAN)[[3]] = dimnames(TAB.RES.MED)[[3]] = 
  paste0("n = ", n.grid)
dimnames(TAB.RES.MEAN)[[4]] = dimnames(TAB.RES.MED)[[4]] = 
  paste0("nsr = ", Nsr.grid)
#
dimnames(RESULTS)[[2]] = dimnames(TAB.RES.MEAN)[[2]]
# save outputs into files
save(TAB.RES.MEAN, file = 
       pathmean<-paste0("data/RESULTS.MEAN.",sim.name,".RData"))
save(TAB.RES.MED, file = paste0("data/RESULTS.MED.",sim.name,".RData"))

##################
#
# Figures
#
##################

(load(pathmean))

all.equal(TAB.RES.MEAN[,1,,],TAB.RES.MEAN[,5,,])
all.equal(TAB.RES.MEAN[,3,,],TAB.RES.MEAN[,7,,])

dimnames(TAB.RES.MEAN)
indices = matrix(c(1,3,2,4,6,8),nrow=2)
fignames = c("reg_bandwidth_choice_rmsep",
             "deriv_bandwidth_choice_rmsep",
             "deriv_bandwidth_choice_hreg_rmsep")
locat = matrix(c(.71,.6625, 0.42-.018, 0.387-.018, 0.3544+.041, 0.3306+.041),
               nrow=2)
for(kk in 1:ncol(indices)){
  pdf(file = paste0("graphics/",fignames[kk],".pdf"))
  RES = TAB.RES.MEAN[1,indices[,kk],,]
  RES.sd = TAB.RES.MEAN[2,indices[,kk],,]
  xpoints = n.grid
  nbxpoints = length(xpoints)
  plot_colors = c("black", "gray") # c("blue","forestgreen")
  yrange = c(0,.74) # range(cbind(RES + RES.sd, RES - RES.sd))
  ylab = expression(ORMSEP[reg])
  if(kk>1){
    yrange = c(0,max(TAB.RES.MEAN[1,indices[,2],,])+0.05) # range(cbind(
      # TAB.RES.MEAN[1,indices[,2],,]+TAB.RES.MEAN[2,indices[,2],,],
      # TAB.RES.MEAN[1,indices[,2],,]-TAB.RES.MEAN[2,indices[,2],,]))
    ylab = expression(ORMSEP[deriv])
  }
  par(mfrow = c(1, 1), cex.axis = 1.5, mar = c(5,5,1,1), cex.lab = 1.5)
  plot(xpoints, RES[1,,1], ylab = ylab, ylim = yrange, xaxt = "n", 
       xlab = "learning sample size", type="n", bty = "L", xlim = c(n.min, 625))
  axis(side = 1, at = xpoints, labels = xpoints)
  for(k in 1:(dim(RES)[3])){
    MED = t(RES[,,k])
    STD = t(RES.sd[,,k])
    lines(xpoints, MED[, 1], type = "b", col = plot_colors[1], pch = 0, 
          lwd = 3, lty = 2)
    lines(xpoints, MED[, 2], type = "b", col = plot_colors[2], pch = 2, 
          lwd = 3, lty = 4)
    # add +/- sd
    # for(jj in c(1,2)){
    #   for(ii in 1:nbxpoints){
    #     segments(xpoints[ii] - 3, MED[ii, jj] - STD[ii, jj], xpoints[ii] + 3, 
    #                 MED[ii, jj] - STD[ii, jj], col = plot_colors[jj], lwd = 1)
    #     # segments(xpoints[ii]+(jj-1.5)*4, MED[ii,jj]-STD[ii,jj], 
    #     #             xpoints[ii]+(jj-1.5)*4, MED[ii,jj]+STD[ii,jj], 
    #     #             col = plot_colors[jj], lwd = 1)
    #     segments(xpoints[ii] - 3, MED[ii, jj] + STD[ii, jj], xpoints[ii] + 3, 
    #                 MED[ii, jj] + STD[ii, jj], col = plot_colors[jj], lwd = 1)
    #   }
    # }
  text(max(xpoints) + 10, MED[nbxpoints, 1]+0.005*(k-2), 
       labels = paste0("nsr = ",Nsr.grid[k]*100,"%"), cex = 1.50, pos = 4)
  }
  legend("topright", legend = c(expression(AIC[C]), "CV"), text.col = 
           plot_colors[c(2,1)], col = plot_colors[c(2,1)], lty = c(4,2), 
         lwd = 3, cex = 2, seg.len = 3, bty = "n")
  points(rep(422, 2), locat[,kk], pch = c(2, 0), col = plot_colors[c(2,1)], 
         lwd = 3)
  dev.off()
}

###################################
#
# Single run of one simulation with figures
#
###################################

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

pdf(file = paste0("graphics/reg_scatterplot_n",fnlearn[1],"_nsr04.pdf"))
par(mfrow = c(1, 1), mar = c(5, 5.35, 1, 1), cex.lab = 2, cex.axis = 1.5)
plot(fRegression[[1]][-(1:(fnlearn[1]))], res[[1]]$Predicted.responses, 
     xlab = expression(m(X[i])), ylab = expression(hat(m)(X[i])), 
     xlim = xbounds, ylim = ybounds)
dev.off()

pdf(file = paste0("graphics/reg_scatterplot_n",fnlearn[2],"_nsr005.pdf"))
par(mfrow = c(1, 1), mar = c(5, 5.35, 1, 1), cex.lab = 2, cex.axis = 1.5)
plot(fRegression[[2]][-(1:(fnlearn[2]))], res[[2]]$Predicted.responses, 
     xlab = expression(m(X[i])), ylab = expression(hat(m)(X[i])), 
     xlim = xbounds, ylim = ybounds)
dev.off()

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

###############################
#
# Comparison of the estimated bandwidth with with oracle bandwidth
#
###############################

B = 100
h.oracle = array(dim=c(length(n.grid), length(Nsr.grid), B))
h.reg = array(dim=c(length(n.grid), length(Nsr.grid), B))
h.deriv = array(dim=c(length(n.grid), length(Nsr.grid), B))
count.n = 0
for(nlearn in n.grid){
  count.n = count.n + 1
  ntest = 2
  sample.size = nlearn + ntest
  count.nsr = 0
  for(nsr in Nsr.grid){
    count.nsr = count.nsr + 1
    cat("n = ", nlearn, " and nsr = ", nsr, "\n", sep ="")
    RESULTS = foreach(b = 1:B, .combine='rbind') %dopar% {
    ################################
    # Simulate all
    ################################
    dat = generate(nlearn, ntest, Jmodel, nsr, model="M1")
    LEARN = dat$LEARN
    PRED = dat$PRED
    Responses = dat$Responses
    oracle.Responses = dat$oracle.Responses
    oracle.DERIV = dat$oracle.DERIV
    DERIV = dat$DERIV  # true derivative functions
    gridsize = ncol(LEARN)
    ##############################################
    # Functional local linear estimator
    ##############################################
    h.ora = oracle.bw(Responses, LEARN, DERIV,BASIS=BASIS)
    res = fllr(Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30, 
               percent = 0.95, BASIS = BASIS, J.est = FALSE)
    return(c(h.ora, res$knn.opt$reg, res$knn.opt$deriv))
    }
    h.oracle[count.n,count.nsr,] = RESULTS[,1]
    h.reg[count.n,count.nsr,] = RESULTS[,2]
    h.deriv[count.n,count.nsr,] = RESULTS[,3]
  }
}

save(h.oracle,h.reg,h.deriv,file=filename<-"data/RESULTS.BWD.RData")

(load(filename))

MEAN.h.oracle = apply(h.oracle,c(1,2),mean)
SD.h.oracle = apply(h.oracle,c(1,2),sd)
MEAN.h.reg = apply(h.reg,c(1,2),mean)
SD.h.reg = apply(h.reg,c(1,2),sd)
MEAN.h.deriv = apply(h.deriv,c(1,2),mean)
SD.h.deriv = apply(h.deriv,c(1,2),sd)

pdf(file = "graphics/deriv_bandwidth_choice.pdf")
par(mfrow = c(1, 1), mar = c(5, 5, 1, 1), cex.lab = 1.5, cex.axis = 1.5)
x = MEAN.h.oracle
y = MEAN.h.deriv
Bounds.x = range(x)
Bounds.y = range(c(MEAN.h.deriv+SD.h.deriv, MEAN.h.reg-SD.h.reg))
symbols(x, y, circles = rep(0.6, length(Nsr.grid)*length(n.grid)), bg = "black",
        inches = FALSE,
        xlab = "oracle bandwidth for the functional derivative",
        ylab = "selected bandwidth for the functional derivative",
        xlim = Bounds.x, ylim = Bounds.y)
sd = SD.h.deriv
segments(x, y - sd, x, y + sd)
width = 0.5
segments(x - width, y - sd, x + width, y - sd)
segments(x - width, y + sd, x + width, y + sd)
abline(0, 1, lty = 2)
#
y = MEAN.h.reg
symbols(x, y, squares = rep(1, length(Nsr.grid)*length(n.grid)), inches = FALSE,
        bg = "gray", xlab = "", ylab ="", xlim = Bounds.x, ylim = Bounds.y, 
        add = T)
sd = SD.h.reg
segments(x, y - sd, x, y + sd, col = "gray")
width = 0.5
segments(x - width, y - sd, x + width, y - sd)
segments(x - width, y + sd, x + width, y + sd)
dev.off()

RES = array(dim=c(2,3,dim(MEAN.h.oracle)))
RES[1,1,,] = MEAN.h.reg
RES[2,1,,] = SD.h.reg
RES[1,2,,] = MEAN.h.deriv
RES[2,2,,] = SD.h.deriv
RES[1,3,,] = MEAN.h.oracle
RES[2,3,,] = SD.h.oracle
dimnames(RES)[[1]] = c("mean", "sd")
dimnames(RES)[[2]] = c("reg", "deriv", "oracle")
dimnames(RES)[[3]] = paste0("n = ", n.grid)
dimnames(RES)[[4]] = paste0("nsr = ", Nsr.grid)

Tab = array(dim = c(3,length(n.grid), length(Nsr.grid)))
dimnames(Tab)[[1]] = c("$h_{reg}$","$h_{deriv}$","$h_{deriv}^{oracle}$")
dimnames(Tab)[[2]] = n.grid
dimnames(Tab)[[3]] = Nsr.grid
for(i in 1:(dim(Tab)[1])) for(j in 1:(dim(Tab)[2])) for(k in 1:(dim(Tab)[3])){
  Tab[i,j,k] = paste(sprintf("%.3f",RES[1,i,j,k])," \\,(",sprintf("%.3f",
                                                       RES[2,i,j,k]),")",sep="")
}
labl = "tab:der_bandwidth_choice"
capt = "Estimated bandwidths."
fTab = ftable(Tab,row.vars=c(2,3),col.vars=1)
xfTab = xtableFtable(fTab, method="compact", label=labl,
                     caption = capt,
                     align=c("c","c","c","c","c","c"), lsep="")
xtbl = print(xfTab,hline.after=c(0,1,1+3*(1:length(n.grid))), 
             caption.placement = "top", table.placement="htpb",
             sanitize.text.function = identity)
write(xtbl,file=filename<-paste("tables/Tab-BW.tex",sep=""))
Tab = readLines(filename)
Tab[7] = sub("&", "$n$ & $nsr$", Tab[7])
write(Tab,filename)