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
TAB.RES.COR = array(NA, c(nbsamples, 2, length(n.grid), length(Nsr.grid)))
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
      ##############################################
      # Functional local linear estimator
      ##############################################
      ### cv
      res = fllr(Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30,
                 percent = 0.95, BASIS = BASIS, J.est = FALSE)
      cor.reg.oracle.cv = cor(oracle.Responses,res$Predicted.responses)
      ### aicc
      res = fllr(Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30,
                 percent = 0.95, BASIS = BASIS, method="aicc", J.est = FALSE)
      cor.reg.oracle.aicc = cor(oracle.Responses,res$Predicted.responses)
      ###
      return(c(cor.reg.oracle.cv, cor.reg.oracle.aicc))
    }
    TAB.RES.COR[, , count.n, count.nsr] = RESULTS
  }
}
#########################################
#
# END LOOP
#
#########################################
dimnames(TAB.RES.COR)[[1]] = paste("sample",1:nbsamples)
dimnames(TAB.RES.COR)[[2]] = c("cor cv", "cor aicc")
dimnames(TAB.RES.COR)[[3]] = paste0("n = ", n.grid)
dimnames(TAB.RES.COR)[[4]] = paste0("nsr = ", Nsr.grid)

save(file="Results_correlation.RData",TAB.RES.COR)

### Figure with correlation boxplots

(load(file="Results_correlation.RData"))

pdf(file = paste0("correlation.pdf"))
par(mfrow = c(1, 1), cex.axis = 1.5, mar = c(5,5,1,1), cex.lab = 1.5)
boxplot(TAB.RES.COR[,2,,1],ylim=c(.9,1),type="n", xaxt = "n",
        xlab = "learning sample size", 
        ylab = "", axes=FALSE)
# title(ylab = expression(cor(hat(m)(X[i]),m(X[i]))),line=2.5)
title(ylab="Pearson's correlation")
axis(side = 1, at = c(-5,1,3,5,7,9), labels = c("",100,200,300,400,500))
axis(side = 1, at = 1:9, labels= rep(NA,9))
axis(side = 2, at = c(-5,.9,.95,1))
dev.off()

#####
##### computation times
#####

n.grid = seq(100, 500, by=50) # sample sizes
Timing = matrix(0,nrow=length(n.grid),ncol=5)
count.n = 0
for(nlearn in n.grid){
  count.n = count.n + 1 
  reps = 10 # number of replicates
  ################################
  # Simulate all
  ################################
  dat = generate(nlearn, ntest = 500, Jmodel = 4, nsr = 0.4, model="M1")
  LEARN = dat$LEARN
  PRED = dat$PRED
  Responses = dat$Responses
  oracle.Responses = dat$oracle.Responses
  oracle.DERIV = dat$oracle.DERIV
  gridsize = ncol(LEARN)
  BASIS = dat$BASIS

  ##############################################
  # Linear model (L)
  ##############################################
  Timing[count.n,1] = 
    system.time({replicate(reps,LinearModel(Responses,LEARN,PRED))})["elapsed"]
  ##############################################
  # Nadaraya-Watson (LC)
  ##############################################
  Timing[count.n,2] = system.time(replicate(reps,{
    DIST = L2metric(LEARN,LEARN)/sqrt(gridsize-1)
    DIST.PRED = L2metric(PRED,LEARN)/sqrt(gridsize-1)
    res.LC = KernelSmoother(Responses,DIST,DIST.PRED, fullgrid=TRUE)
  }))["elapsed"]
  ##############################################
  # Functional local linear estimator (LL)
  ##############################################
  Timing[count.n,3] = system.time(
    replicate(reps,{fllr(Responses, LEARN, PRED)}))["elapsed"]
  ##############################################
  # Mueller-Yao method (FAM)
  ##############################################
  Timing[count.n,4] = system.time({
    replicate(reps, FAM(Responses,LEARN, PRED, fullgrid=TRUE))
    })["elapsed"]
  ##############################################
  # Gaussian process regression
  ##############################################
  Timing[count.n,5] = system.time(replicate(reps,{GPR(Responses,LEARN, PRED)}))["elapsed"]
}

dimnames(Timing)[[1]] = n.grid
dimnames(Timing)[[2]] = c("L","LC","LL","FAM","GPR")
RES = Timing/reps

labl = paste0("tab:timing")
capt = paste0("Average running times (in s.) of the competing methods, 
              based on 10 replications of the experiment with learning sample size $n$.")
write(print(xtable(RES,align=c("c|","c","c","c","c","c"),label=labl, 
                   caption = capt, digits=3),hline.after=c(0), 
            caption.placement = "top", table.placement="htpb",
            sanitize.text.function = identity),file=
        filename<-paste0("Tab_timing.tex"))
Tab = readLines(filename)
Tab[9] = paste0("$n$", Tab[9])
write(Tab,file=filename)
