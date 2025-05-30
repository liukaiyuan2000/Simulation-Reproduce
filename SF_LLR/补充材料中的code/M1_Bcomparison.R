################################################################################
# 
# Ferraty, F. and Nagy, S. (2021). Scalar-on-function local linear regression 
#                                  and beyond.
#
#  Model (M1)
#  J=4 and the Fourier basis for building the regression model
#  Fourier basis is used to predict the functions
#  data-driven choice of 
#		* two bandwidths: - one for the regression operator and 
#                     - one for its functional derivatives
#
################################################################################

for(Bsetup in 1:2){
sim.name = paste0("M1_B",Bsetup)  # identifier of the simulation study
if(Bsetup == 1){
  nsr = 0.4
  nlearn = n.min
}
if(Bsetup == 2){
  nsr = .05
  nlearn = n.max
}

################################################################################
# Getting the distribution of kNN bandwidths and quality criteria with respect 
# to the sample size / noise-to-signal ratio / bandwidth selection method
################################################################################
set.seed(1)
nbsamples = NBsamples       # number of independent runs
Jmodel = 4                  # number of basis elements used 
                            # for building functional predictors
#########################################
#
# START LOOP
#
#########################################

sample.size = nlearn + ntest
B.grid = c(50,100,500,1000)
lB = length(B.grid)
TAB.RES.MED = array(NA, c(2, 6, lB))
TAB.RES.MEAN = array(NA, c(2, 6, lB))
##
cat("Bsetup = ", Bsetup,"\n", sep ="")
acomb <- function(...) abind(..., along=3)
RESULTS = foreach(kk = 1:nbsamples, .combine='acomb') %dopar% {
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
    count.B = 0
    msep.reg.oracle = msep.deriv.oracle = knn.opt.reg = knn.opt.deriv = 
      min.crit.reg = dim.opt = rep(NA,lB)
    for(B in B.grid){
      count.B = count.B + 1
      res = fllr(Responses, LEARN, PRED, nboot = B.grid[count.B], 
                 kNN.grid.length = 30,
                 percent = 0.95, BASIS = BASIS, J.est = FALSE, 
                 fullgrid = TRUE)
      msep.reg.oracle[count.B] = 
        mean((oracle.Responses - res$Predicted.responses)^2)/denom.reg
      msep.deriv.oracle[count.B] = 
        mean((oracle.DERIV - res$PREDICTED.DERIV)^2)/denom.deriv
      knn.opt.reg[count.B] = res$knn.opt$reg
      knn.opt.deriv[count.B] = res$knn.opt$deriv
      min.crit.reg[count.B] = res$min.crit.reg
      dim.opt[count.B] = res$dim.opt
    }
    ###
    return(cbind(knn.opt.reg, knn.opt.deriv, min.crit.reg, 
             msep.reg.oracle, msep.deriv.oracle, dim.opt
    ))
}
    for(count.B in 1:lB){
      TAB.RES.MED[, , count.B] = t(cbind(
        apply(RESULTS[count.B,,], 1 , median), 
        apply(RESULTS[count.B,,], 1, mad)))
      TAB.RES.MEAN[, , count.B] = t(cbind(
        apply(RESULTS[count.B,,], 1, mean),
        apply(RESULTS[count.B,,], 1, sd)))
    }
#########################################
#
# END LOOP
#
#########################################
dimnames(TAB.RES.MEAN)[[1]] = c("mean", "sd")
dimnames(TAB.RES.MED)[[1]] = c("median", "mad")
dimnames(TAB.RES.MEAN)[[2]] = dimnames(TAB.RES.MED)[[2]] = c(
  "knn.reg", "knn.deriv", "mse reg", 
  "msep reg oracle", "msep deriv oracle", "optimal dimension")
dimnames(TAB.RES.MEAN)[[3]] = dimnames(TAB.RES.MED)[[3]] = paste0("B = ",B.grid)
#
dimnames(RESULTS)[[2]] = dimnames(TAB.RES.MEAN)[[2]]
# save outputs into files
save(TAB.RES.MEAN, file = 
       pathmean<-paste0("data/RESULTS.MEAN.",sim.name,".RData"))
save(TAB.RES.MED, file = paste0("data/RESULTS.MED.",sim.name,".RData"))
}

####################
#
# Table of results
#
####################

for(Bsetup in 1:2){
  sim.name = paste0("M1_B",Bsetup)  # identifier of the simulation study
  if(Bsetup == 1){
    nsr = 0.4
    nlearn = n.min
  }
  if(Bsetup == 2){
    nsr = .05
    nlearn = n.max
  }
B.grid = c(50,100,500,1000)
(load(pathmean<-paste0("data/RESULTS.MEAN.",sim.name,".RData")))
RES = t(TAB.RES.MEAN[1,c(2,5),])
RES.sd = t(TAB.RES.MEAN[2,c(2,5),])
Tab = matrix(NA,nrow(RES),ncol(RES))
rownames(Tab) = B.grid
colnames(Tab) = c("$h_{deriv}$", "$ORMSEP_{deriv}$")
for(i in 1:(dim(Tab)[1])) for(j in 1:(dim(Tab)[2])){
  Tab[i,j] = paste(sprintf("%.3f",RES[i,j])," \\,(",
                   sprintf("%.3f",RES.sd[i,j]),")",sep="")
}
if(Bsetup == 1) labl = "TableS1(a)"
if(Bsetup == 2) labl = "TableS1(b)"
capt = paste0("Stability of the bootstrap bandwidth selection, $n = ",
              nlearn ,"$ and $nsr = ", nsr, "$.")
write(print(xtable(Tab,align=c("c|","c","c"),label=labl, caption = capt),
            hline.after=c(0), caption.placement = "top", table.placement="htpb",
            sanitize.text.function = identity),file=
        filename<-paste0("tables/TabS1_",Bsetup,".tex"))
Tab = readLines(filename)
Tab[8] = paste0("$B$", Tab[8])
write(Tab,filename)
}

## compile the two obtained tables into a single one
mainfilename = paste0("tables/TabS1.tex")
filename1  = paste0("tables/TabS1_1.tex")
filename2  = paste0("tables/TabS1_2.tex")
Tab = readLines(mainfilename)
Tab1 = readLines(filename1) 
Tab2 = readLines(filename2) 
Tab[13:16] = Tab1[10:13]
Tab[23:26] = Tab2[10:13]
write(Tab,mainfilename)
