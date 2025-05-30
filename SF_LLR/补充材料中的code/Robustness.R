################################################################################
# 
# Ferraty, F. and Nagy, S. (2021). Scalar-on-function local linear regression 
#                                  and beyond.
#
#  Robustness study
#  J=2,3,4 and the Fourier basis for building the regression model
#  Fourier basis is used to predict the functions
#  data-driven choice of 
#		* two bandwidths: - one for the regression operator and 
#                     - one for its functional derivatives
#
################################################################################
# Getting the distribution of kNN bandwidths and quality criteria with respect 
# to the sample size / noise-to-signal ratio / bandwidth selection method
################################################################################
set.seed(1)
nbsamples = NBsamples       # number of independent runs
Jmodel = 4                  # maximum number of basis elements used 
                            # for building functional predictors
#########################################
#
# START LOOP
#
#########################################

n.grid = seq(n.min, n.max, by=50) # sample sizes
Nsr.grid = c(0.05, 0.4)
rho.grid = c(0.05,0.1,0.2,0.4)
J.grid = c(2,3,4)
ncolTAB = 6
TAB.RES.MED = array(NA, c(2, ncolTAB, length(n.grid), length(Nsr.grid), 
                          length(rho.grid), length(J.grid)))
TAB.RES.MEAN = array(NA, c(2, ncolTAB, length(n.grid), length(Nsr.grid), 
                           length(rho.grid), length(J.grid)))
##
count.n = 0
for(nlearn in n.grid){
  count.n = count.n + 1
  sample.size = nlearn + ntest
  count.nsr = 0
  for(nsr in Nsr.grid){
    count.nsr = count.nsr + 1
    count.rho = 0
    for(rho in rho.grid){
      count.rho = count.rho + 1
      count.J = 0
      for(J in J.grid){
        count.J = count.J + 1
        cat("n = ", nlearn, ", nsr = ", nsr,", rho = ", rho, ", J = ",
            J,"\n", sep ="")
        RESULTS = foreach(kk = 1:nbsamples, .combine='rbind') %dopar% {
        ################################
        # Simulate all
        ################################
        dat = generate(nlearn, ntest, Jmodel=J, nsr=nsr, rho=rho, model="M2")
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
        res = fllr(Responses, LEARN, PRED, nboot = 100, kNN.grid.length = 30, 
                   percent = 0.95)
        msep.reg.oracle = mean((oracle.Responses - res$Predicted.responses)^2)/
          denom.reg
        msep.deriv.oracle = mean((oracle.DERIV - res$PREDICTED.DERIV)^2)/
          denom.deriv
        ##
        return(c(res$knn.opt$reg, res$knn.opt$deriv, res$min.crit.reg, 
                 msep.reg.oracle, msep.deriv.oracle, res$dim.opt==J
               ))
        }
      TAB.RES.MED[, , count.n, count.nsr, count.rho, count.J] = 
        t(cbind(apply(RESULTS, 2, median), apply(RESULTS, 2, mad)))
      TAB.RES.MEAN[, , count.n, count.nsr, count.rho, count.J] = 
        t(cbind(apply(RESULTS, 2, mean), apply(RESULTS, 2, sd)))
      }
    }
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
  "knn.reg", "knn.deriv", 
  "mse reg", "msep reg oracle", "msep deriv oracle", "optimal dimension")
dimnames(TAB.RES.MEAN)[[3]] = dimnames(TAB.RES.MED)[[3]] = 
  paste0("$n = ", n.grid,"$")
dimnames(TAB.RES.MEAN)[[4]] = dimnames(TAB.RES.MED)[[4]] = 
  paste0("nsr = ", Nsr.grid)
dimnames(TAB.RES.MEAN)[[5]] = dimnames(TAB.RES.MED)[[5]] = 
  paste0("$\\rho = ", rho.grid,"$")
dimnames(TAB.RES.MEAN)[[6]] = dimnames(TAB.RES.MED)[[6]] = 
  paste0("$J = ", J.grid,"$")
#
dimnames(RESULTS)[[2]] = dimnames(TAB.RES.MEAN)[[2]]
# save outputs into files
save(TAB.RES.MEAN, file = 
       pathmean<-paste0("data/RESULTS.MEAN.",sim.name,".RData"))
save(TAB.RES.MED, file = paste0("data/RESULTS.MED.",sim.name,".RData"))

####################
#
# Tables of results
#
####################

(load(pathmean))
TAB = TAB.RES.MEAN

library(xtable)
# Table 1: number of correctly selected dimension
# Table 2: ORMSEP_reg
# Table 3: ORMSEP_deriv
# Table 4: ORMSEP both as ftable
Nsr.grid.f = c("005","04")
for(count.nsr in 1:2) for(count.J in 1:length(J.grid)){
  # number of correct dimensions
  RES = TAB[1,6,,count.nsr,,count.J]*NBsamples
  rownames(RES) = n.grid
  labl = paste0("tab:dimension_choice_J",J.grid[count.J],"_nsr",
                Nsr.grid.f[count.nsr])
  capt = paste0("Number of times, out of $",NBsamples,"$, 
                that the dimension is correctly selected with $J=",
                J.grid[count.J],"$, $nsr=",Nsr.grid[count.nsr],"$.")
  write(print(xtable(RES,align=c("c|","c","c","c","c"),label=labl, caption = 
                       capt, digits=0),hline.after=c(0), caption.placement = 
                "top", table.placement="htpb",
              sanitize.text.function = identity),file=
          filename<-paste0("tables/Tab1_",count.nsr,count.J,".tex"))
  Tab = readLines(filename)
  Tab[8] = paste0("$n$", Tab[8])
  write(Tab,file=filename)
  # ORMSEP_reg
  RES = TAB[1,4,,count.nsr,,count.J]
  RES.sd = TAB[2,4,,count.nsr,,count.J]
  Tab = matrix(NA,nrow(RES),ncol(RES))
  for(i in 1:nrow(RES)) for(j in 1:ncol(RES)){
    Tab[i,j] = paste(sprintf("%.3f",RES[i,j])," \\,(",
                     sprintf("%.3f",RES.sd[i,j]),")",sep="")
  }
  rownames(Tab) = n.grid
  colnames(Tab) = colnames(RES)
  labl = paste0("tab:ormsep_reg_J",J.grid[count.J],"_nsr",Nsr.grid.f[count.nsr])
  capt = paste0("Average and standard deviation (in brackets) of 
                $ORMSEP_{reg}$ with $J=",J.grid[count.J],"$, $nsr=",
                Nsr.grid[count.nsr],"$.")
  write(print(xtable(Tab,align=c("c|","c","c","c","c"),label=labl, 
                     caption = capt),hline.after=c(0), 
              caption.placement = "top", table.placement="htpb",
              sanitize.text.function = identity),file=
          filename<-paste0("tables/Tab2_",count.nsr,count.J,".tex"))
  Tab = readLines(filename)
  Tab[8] = paste0("$n$", Tab[8])
  write(Tab,file=filename)
  # ORMSEP_deriv
  RES = TAB[1,5,,count.nsr,,count.J]
  RES.sd = TAB[2,5,,count.nsr,,count.J]
  Tab = matrix(NA,nrow(RES),ncol(RES))
  for(i in 1:nrow(RES)) for(j in 1:ncol(RES)){
    Tab[i,j] = paste(sprintf("%.3f",RES[i,j])," \\,(",
                     sprintf("%.3f",RES.sd[i,j]),")",sep="")
  }
  rownames(Tab) = n.grid
  colnames(Tab) = colnames(RES)
  labl = paste0("tab:ormsep_deriv_J",J.grid[count.J],"_nsr",
                Nsr.grid.f[count.nsr])
  capt = paste0("Average and standard deviation (in brackets) of 
                $ORMSEP_{deriv}$ with $J=",J.grid[count.J],"$, $nsr=",
                Nsr.grid[count.nsr],"$.")
  write(print(xtable(Tab,align=c("c|","c","c","c","c"),label=labl, 
                     caption = capt),hline.after=c(0), 
              caption.placement = "top", table.placement="htpb",
              sanitize.text.function = identity),file=
          filename<-paste0("tables/Tab3_",count.nsr,count.J,".tex"))
  Tab = readLines(filename)
  Tab[8] = paste0("$n$", Tab[8])
  write(Tab,file=filename)
  # ORMSEP both
  RES = TAB[1,c(4,5),,count.nsr,,count.J]
  RES.sd = TAB[2,c(4,5),,count.nsr,,count.J]
  Tab = array(dim=dim(RES))
  for(i in 1:(dim(RES)[1])) for(j in 1:(dim(RES)[2])) for(k in 1:(dim(RES)[3])){
    Tab[i,j,k] = paste(sprintf("%.3f",RES[i,j,k])," \\,(",
                       sprintf("%.3f",RES.sd[i,j,k]),")",sep="")
  }
  dimnames(Tab)[[1]] = c(
    paste0("\\parbox[t]{2mm}{\\multirow{",length(n.grid),
           "}{*}{\\rotatebox[origin=c]{90}{$ORMSEP_{reg}$}}}"),
    paste0("\\parbox[t]{2mm}{\\multirow{",length(n.grid),
           "}{*}{\\rotatebox[origin=c]{90}{$ORMSEP_{deriv}$}}}"))
  dimnames(Tab)[[2]] = dimnames(RES)[[2]]
  dimnames(Tab)[[3]] = dimnames(RES)[[3]]
  labl = paste0("tab:ormsep_J",J.grid[count.J],"_nsr",Nsr.grid.f[count.nsr])
  capt = paste0("Average and standard deviation (in brackets) of $ORMSEP$ with 
                $J=",J.grid[count.J],"$, $nsr=",Nsr.grid[count.nsr],"$.")
  xtbl = xtableFtable(ftable(Tab,row.vars=c(1,2),col.vars=3),method="compact",
               label=labl,
               caption =capt,
               align=c("c","c","c|","c","c","c","c"),
               lsep="")
  xtblp = print(xtbl,
                hline.after=c(1,length(n.grid)+1),
                caption.placement = "top", table.placement="htpb")
  write(xtblp, file = filename<-paste0("tables/Tab4_",count.nsr,count.J,".tex"))
  Tab = readLines(filename)
  Tab[6] = sub("&", "& $n$ ", Tab[6])
  Tab = gsub("n = ","",Tab)
  write(Tab,file=filename)
}

###############
#
# Table with computation times
#
###############

count.nsr = 1
count.J = length(J.grid)
n.grid = seq(n.min, n.max, by=50) # sample sizes
nsr = 0.05
rho = 0.05
J = 4
Timing = rep(NA,length(n.grid))
count.n = 0
for(nlearn in n.grid){
  count.n = count.n + 1
  sample.size = nlearn + ntest
  ################################
  # Simulate all
  ################################
  dat = generate(nlearn, ntest, Jmodel=J, nsr=nsr, rho=rho, model="M2")
  LEARN = dat$LEARN
  PRED = dat$PRED
  Responses = dat$Responses
  ##############################################
  # Functional local linear estimator
  ##############################################
  Timing[count.n] = system.time({res = fllr(Responses, LEARN, PRED, 
                 nboot = 100, kNN.grid.length = 30, percent = 0.95)})["elapsed"]
}

# number of times correct dimension was selected, and timings for the main paper
RES = TAB[1,6,,count.nsr,,count.J]*NBsamples
rownames(RES) = n.grid
RES = cbind(RES,Timing)
labl = paste0("tab:dimension_choice")
capt = paste0("Number of times, out of $",NBsamples,"$, that the dimension is 
              correctly selected, and the running times (in s.).")
dgts = matrix(0,nrow=nrow(RES),ncol=ncol(RES)+1)
dgts[,ncol(RES)+1] = 2
write(print(xtable(RES,align=c("c|","c","c","c","c", "|c"),label=labl, 
                   caption = capt, digits=dgts),hline.after=c(0), 
            caption.placement = "top", table.placement="htpb",
            sanitize.text.function = identity),file=
        filename<-paste0("Tab1_main.tex"))
Tab = readLines(filename)
Tab[8] = paste0("$n$", Tab[8])
write(Tab,file=filename)