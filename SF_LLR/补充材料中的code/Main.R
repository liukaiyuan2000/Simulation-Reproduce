################################################################################
# 
# Ferraty, F. and Nagy, S. (2021). Scalar-on-function local linear regression 
#                                  and beyond.
# 
# Main file for the simulation studies, from which tasks are distributed
#
################################################################################
rm(list=ls())
library(xtable)
library(fllr)    # package available in on-line Supplementary Material
library(abind)   # combining multidimensional arrays along dimensions
NBsamples = 100  # overall number of runs in all studies
n.min = 100
ntest = 500      # learning sample size

#######################################
# Setting parallelization features
#######################################
library(doSNOW)
### set 10 parallel jobs; you can adapt this number of job
### according to the power of your computer; if you have no
### idea about its capacity, start with 1 job, and try again
### with 2 jobs, and so on (with a too large number of parallel
### jobs the ram of your computer may be saturated and then
### the computation may slow down dramatically!)
nbofjobs = 10
cl = makeCluster(rep("localhost", nbofjobs), type = "SOCK")
registerDoSNOW(cl)
clusterEvalQ(cl, library(doSNOW))
clusterEvalQ(cl, library(parallelDist))
clusterEvalQ(cl, library(fllr))
# this configuration may be canceled at any time with the R command 
# "stopCluster(cl)"

Jmodel = 4    # number of basis elements used for building functional predictors
tm1 = proc.time()
#
n.max = 100
sim.name = "M2_100_J4"
source("M2.R")
#
tm2 = proc.time()
n.max = 500
sim.name = "M2_500_J4"
source("M2.R")

Jmodel = 15   # number of basis elements used for building functional predictors
tm3 = proc.time()
#
n.max = 100
sim.name = "M2_100_J15"
source("M2.R")
#
tm4 = proc.time()
n.max = 500
sim.name = "M2_500_J15"
source("M2.R")

n.max = 500
#
tm5 = proc.time()
sim.name = "Robustness"
source("Robustness.R")   
# complete study with (M1*, M1 with rho) - not in the paper, 
# output are:
#   Tab1_ij.tex; Tab2_ij.tex; Tab3_ij.tex 
#   for i=1,2 (nsr level) and j=1,2,3 (J-1 dimension)

tm6 = proc.time()
source("M1_Bcomparison.R")  # Table S1
# output are:
#   TabS1_1; TabS1_2
tm7 = proc.time()
sim.name = "M1_Figures"
source("M1_Figures.R")   # Figures with comparison of CV and AICC
# further output:
#   Tab-BW.tex

tm8 = proc.time()
source("growth_dataset.R") # Figures for the real data example
stopCluster(cl)
tm8-tm5