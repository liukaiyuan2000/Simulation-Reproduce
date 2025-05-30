################################################################################
# 
# Ferraty, F. and Nagy, S. (2021). Scalar-on-function local linear regression 
#                                  and beyond.
#
#  Model (M2): Main simulation study
#  J=4 or J=15 and the Fourier basis for building the regression model
#  eigenbasis is used to predict the functions
#  data-driven choice of 
#		* two bandwidths: - one for the regression operator and 
#                     - one for its functional derivatives
#   * dimension of the eigenspace to be used
#
################################################################################

set.seed(1)
nlearn = n.max                           # testing sample size
sample.size = nlearn + ntest             # whole sample size
nbsamples = NBsamples                    # number independent runs
#########################################
#
# START LOOP
#
#########################################
Nsr.grid = c(0.05, 0.1, 0.2, 0.4)
Coef.grid = c(0, 0.25, 0.5, 0.75, 1) 
# coefficient balancing the linear part with the nonparametric one
# OMSEP versions
ncolTAB = 13+4
TAB.RES.MED = array(NA, c(2, ncolTAB, length(Coef.grid), length(Nsr.grid)))
TAB.RES.MEAN = array(NA, c(2, ncolTAB, length(Coef.grid), length(Nsr.grid)))
# ORMSEP versions
TAB.RES.MED.R = array(NA, c(2, ncolTAB, length(Coef.grid), length(Nsr.grid)))
TAB.RES.MEAN.R = array(NA, c(2, ncolTAB, length(Coef.grid), length(Nsr.grid)))
count.coef = 0
for(aa in Coef.grid){
	count.coef = count.coef + 1
	count.nsr = 0
	for(nsr in Nsr.grid){
		count.nsr = count.nsr + 1
		cat("coef = ", aa, " and nsr = ", nsr, "\n", sep ="")
		RESULTS = foreach(kk = 1:nbsamples, .combine='rbind') %dopar% {
		  ################################
		  # Simulate all
		  ################################
      dat = generate(nlearn, ntest, Jmodel, nsr, aa)
			LEARN = dat$LEARN
			PRED = dat$PRED
			Responses = dat$Responses
			oracle.Responses = dat$oracle.Responses
			oracle.DERIV = dat$oracle.DERIV
			gridsize = ncol(LEARN)
			######################################
			# Variances for denominators of ORMSEP
			######################################
			denom.reg = mean((oracle.Responses-mean(oracle.Responses))^2)
			denom.deriv = mean((t(oracle.DERIV) - colMeans(oracle.DERIV))^2)
			##############################################
			# Functional local linear estimator (LL)
			##############################################
			res = fllr(Responses, LEARN, PRED)
			msep.reg.oracle = mean((oracle.Responses - res$Predicted.responses)^2)
			msep.deriv.oracle = mean((oracle.DERIV - res$PREDICTED.DERIV)^2)
			##############################################
			# Mueller-Yao method (FAM)
			##############################################
			res.FAM = FAM(Responses,LEARN, PRED, fullgrid=TRUE)
			msep.reg.oracle.FAM = mean((oracle.Responses - 
			                             res.FAM$Predicted.responses)^2)
			msep.deriv.oracle.FAM = mean((oracle.DERIV - res.FAM$PREDICTED.DERIV)^2)
			##############################################
			# Nadaraya-Watson (LC)
			##############################################
			DIST = L2metric(LEARN,LEARN)/sqrt(gridsize-1)
			DIST.PRED = L2metric(PRED,LEARN)/sqrt(gridsize-1)
			res.LC = KernelSmoother(Responses,DIST,DIST.PRED, fullgrid=TRUE)
			msep.oracle.LC = mean((oracle.Responses - res.LC$Predicted.responses)^2)
			##############################################
			# Linear model (L)
			##############################################
			res.L = LinearModel(Responses,LEARN,PRED)
			msep.oracle.L = mean((oracle.Responses - res.L$Predicted.responses)^2)
			msep.deriv.oracle.L = mean((oracle.DERIV - res.L$PREDICTED.DERIV)^2)
			##############################################
			# Mueller-Yao method - BIC used for dimension selection (FAMb)
			##############################################
			res.FAMb = FAM(Responses,LEARN, PRED, method="BIC", fullgrid=TRUE)
			msep.reg.oracle.FAMb = mean((oracle.Responses - 
			                              res.FAMb$Predicted.responses)^2)
			msep.deriv.oracle.FAMb = mean((oracle.DERIV - res.FAMb$PREDICTED.DERIV)^2)
			##############################################
			# Gaussian process regression
			##############################################
			res.GPR = GPR(Responses,LEARN, PRED)
			msep.reg.oracle.GPR = mean((oracle.Responses - 
			                              res.GPR$Predicted.responses)^2)
			###
			return(c(res$knn.opt$reg, res$knn.opt$deriv, res$mse,
			         msep.reg.oracle, msep.deriv.oracle, res$dim.opt,
			         msep.oracle.LC,
			         msep.oracle.L, msep.deriv.oracle.L, res.L$dim.opt,
			         msep.reg.oracle.FAM, msep.deriv.oracle.FAM, res.FAM$dim,
			         msep.reg.oracle.FAMb, msep.deriv.oracle.FAMb, res.FAMb$dim,
			         msep.reg.oracle.GPR,
			         denom.reg, denom.deriv
			         ))
		}
	RESULTS0 = RESULTS
	RESULTS = RESULTS[,1:ncolTAB]
  # OMSEP
	TAB.RES.MED[, , count.coef, count.nsr] = t(cbind(apply(RESULTS, 2, median), 
	                                                 apply(RESULTS, 2, mad)))
	TAB.RES.MEAN[, , count.coef, count.nsr] = t(cbind(apply(RESULTS, 2, mean), 
	                                                  apply(RESULTS, 2, sd)))
  # ORMSEP
	ind.reg = c(4,7,8,11,14,17) # indices to be divided by denom.reg
	for(ir in ind.reg) RESULTS0[,ir] = RESULTS0[,ir]/RESULTS0[,ncolTAB+1]
  ind.deriv = c(5,9,12,15)    # indices to be divided by denom.deriv
	for(ir in ind.deriv) RESULTS0[,ir] = RESULTS0[,ir]/RESULTS0[,ncolTAB+2]
	RESULTS = RESULTS0[,1:ncolTAB]
	TAB.RES.MED.R[, , count.coef, count.nsr] = t(cbind(apply(RESULTS, 2, median), 
	                                                   apply(RESULTS, 2, mad)))
	TAB.RES.MEAN.R[, , count.coef, count.nsr] = t(cbind(apply(RESULTS, 2, mean), 
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
  "knn.reg", "knn.deriv", "mse reg", "msep reg oracle", "msep deriv oracle", 
  "optimal dimension",
  "msep oracle LC",
  "msep oracle L", "msep deriv oracle L", "optimal dimension L",
  "msep reg oracle FAM", "msep deriv oracle FAM", "dimension FAM",
  "msep reg oracle FAM_BIC", "msep deriv oracle FAM_BIC", "dimension FAM_BIC",
  "msep reg oracle GPR"
  )
dimnames(TAB.RES.MEAN)[[3]] = dimnames(TAB.RES.MED)[[3]] = 
  paste0("a = ", Coef.grid)
dimnames(TAB.RES.MEAN)[[4]] = dimnames(TAB.RES.MED)[[4]] = 
  paste0("nsr = ", Nsr.grid)
dimnames(TAB.RES.MEAN.R) = dimnames(TAB.RES.MEAN)
dimnames(TAB.RES.MED.R) = dimnames(TAB.RES.MED)
#
dimnames(RESULTS)[[2]] = dimnames(TAB.RES.MEAN)[[2]]
# save outputs into files
save(TAB.RES.MEAN, file = 
       pathmean<-paste0("data\RESULTS.MEAN.",sim.name,".RData"))
save(TAB.RES.MED, file = paste0("data\RESULTS.MED.",sim.name,".RData"))
save(TAB.RES.MEAN.R, file = 
       pathmean.R<-paste0("data\RESULTS.MEAN.R.",sim.name,".RData"))
save(TAB.RES.MED.R, file = paste0("data\RESULTS.MED.R.",sim.name,".RData"))

##################
#
# Tables with results
#
##################

(load(pathmean))
TAB = TAB.RES.MEAN
(load(pathmean.R))
TAB.R = TAB.RES.MEAN.R

library(xtable)
for(k in 1:length(Nsr.grid)){
  RES = TAB[1,ord<-c(8,7,4,11,14,17,9,5,12,15),,k]
  RES.sd = TAB[2,ord,,k]
  #
  Tab = matrix(NA,nrow(RES),ncol(RES))
  for(i in 1:nrow(RES)) for(j in 1:ncol(RES)){
    Tab[i,j] = paste(sprintf("%.3f",RES[i,j])," \\,(",
                     sprintf("%.3f",RES.sd[i,j]),")",sep="")
  }
  dimnames(Tab)[[2]] = paste("$a = ",Coef.grid,"$",sep="")
  dimnames(Tab)[[1]] = c("L", "LC", "LL", "FAM", "FAM2", "GPR",
                         "Lder", "LLder", "FAMder", "FAM2der")
  labl = paste0("tab:",sim.name,"-",k)
  capt = paste("Model (M3) with $nsr = ",Nsr.grid[k],"$ and $\\rho=",
               Nsr.grid[k],"$.",sep="")
  write(print(xtable(Tab,align=c("c|c|","c","c","c","c","c"),label=labl, 
                     caption = capt),hline.after=c(0,6), 
              caption.placement = "top", table.placement="htpb",
            sanitize.text.function = identity),file=
          filename<-paste("tables\Tab-",sim.name,"-",k,".tex",sep=""))
  Tab = readLines(filename)
  for(i in c(8,10:15,17:20)) Tab[i] = paste("& ",Tab[i],sep="")
  Tab[10] = paste(
    "\\parbox[t]{2mm}{\\multirow{6}{*}{\\rotatebox[origin=c]{90}{Reg.}}}",
    Tab[10],sep="")
  Tab[17] = paste(
    "\\parbox[t]{2mm}{\\multirow{4}{*}{\\rotatebox[origin=c]{90}{Deriv.}}}",
    Tab[17],sep="")
  Tab[7] = paste("\\resizebox{\\textwidth}{!}{", Tab[7],sep="")
  Tab[21] = paste(Tab[21], "}",sep="")
  Tab = gsub("der", "", Tab,fixed=TRUE)
  write(Tab,file=filename)
}
###

for(k in 1:length(Nsr.grid)){
  RES = TAB.R[1,ord,,k]
  RES.sd = TAB.R[2,ord,,k]
  #
  Tab = matrix(NA,nrow(RES),ncol(RES))
  for(i in 1:nrow(RES)) for(j in 1:ncol(RES)){
    Tab[i,j] = paste(sprintf("%.3f",RES[i,j])," \\,(",
                     sprintf("%.3f",RES.sd[i,j]),")",sep="")
  }
  # for a = 0 and derivative, put results for OMSEP, in gray
  RES.OMSEP = TAB[1,ord,,k]
  RES.sd.OMSEP = TAB[2,ord,,k]
  for(i in 7:10) Tab[i,1] = paste("\\graytext{",sprintf("%.3f",RES.OMSEP[i,j]),
                          " \\,(",sprintf("%.3f",RES.sd.OMSEP[i,j]),")}",sep="")
  #
  dimnames(Tab)[[2]] = paste("$a = ",Coef.grid,"$",sep="")
  dimnames(Tab)[[1]] = c("L", "LC", "LL", "FAM", "FAM2", "GPR",
                         "Lder", "LLder", "FAMder", "FAM2der")
  labl = paste0("tab:R-",sim.name,"-",k)
  capt = paste("Model (M3) with $nsr = ",Nsr.grid[k],"$ and $\\rho=",
               Nsr.grid[k],"$.",sep="")
  write(print(xtable(Tab,align=c("c|c|","c","c","c","c","c"),label=labl, 
                     caption = capt),hline.after=c(0,6), 
              caption.placement = "top", table.placement="htpb",
              sanitize.text.function = identity),file=
          filename<-paste("tables\Tab.R-",sim.name,"-",k,".tex",sep=""))
  Tab = readLines(filename)
  for(i in c(8,10:15,17:20)) Tab[i] = paste("& ",Tab[i],sep="")
  Tab[10] = paste(
    "\\parbox[t]{2mm}{\\multirow{6}{*}{\\rotatebox[origin=c]{90}{Reg.}}}",
    Tab[10],sep="")
  Tab[17] = paste(
    "\\parbox[t]{2mm}{\\multirow{4}{*}{\\rotatebox[origin=c]{90}{Deriv.}}}",
    Tab[17],sep="")
  Tab[7] = paste("\\resizebox{\\textwidth}{!}{", Tab[7],sep="")
  Tab[21] = paste(Tab[21], "}",sep="")
  Tab = gsub("der", "", Tab,fixed=TRUE)
  write(Tab,file=filename)
}

################
#
# Table with dimensions for TAB.R
#
################

# pathmean.R<-paste0("tables\RESULTS.MEAN.R.",sim.name,".RData")
# (load(pathmean.R))
# TAB.R = TAB.RES.MEAN.R
for(k in 1:length(Nsr.grid)){
  RES = TAB.R[1,ord<-c(10,6,13,16),,k]
  RES.sd = TAB.R[2,ord,,k]
  Tab = matrix(NA,nrow(RES),ncol(RES))
  for(i in 1:nrow(RES)) for(j in 1:ncol(RES)){
    Tab[i,j] = paste(sprintf("%.3f",RES[i,j])," \\,(",
                     sprintf("%.3f",RES.sd[i,j]),")",sep="")
  }
  dimnames(Tab)[[2]] = paste("$a = ",Coef.grid,"$",sep="")
  dimnames(Tab)[[1]] = c("L", "LL", "FAM", "FAM2")
  labl = paste0("tab:dim-",sim.name,"-",k)
  capt = paste("Estimated dimensions: Model (M3) with $nsr = ",
               Nsr.grid[k],"$ and $\\rho=",Nsr.grid[k],"$.",sep="")
  write(print(xtable(Tab,align=c("c|","c","c","c","c","c"),label=labl, 
                     caption = capt),hline.after=c(0), 
              caption.placement = "top", table.placement="htpb",
              sanitize.text.function = identity),file=
          filename<-paste("tables\dim-",sim.name,"-",k,".tex",sep=""))
}