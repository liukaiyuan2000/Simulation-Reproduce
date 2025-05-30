rm(list = ls())
library(kernlab)
library(energy)
library(sm)
source("D:/Desktop/JASAIndependenceTest/8014586/source.R")
source("D:/Desktop/JASAIndependenceTest/code/DGP.R")
Rcpp::sourceCpp("D:/Desktop/JASAIndependenceTest/code/function.cpp")
B = 199
RES = NULL
slice.num.choose = c(2, 5, 10, 23, 46, 115)
length.slice = length(slice.num.choose)
res.slice  = rep(0, length.slice)
res.disco  = rep(0, length.slice)
for(i in 1:length.slice){
  slice.num = slice.num.choose[i]
  cat('Example 4: ', 'Slice number = ', slice.num, '\r', sep = "")
  dat = DGP(example = 4, slice.num = slice.num)
  x = dat$x
  y.d = dat$y.d
  res.slice[i] <- ptest_slm_rcpp(x, y.d, B = B)$pvalue
  res.disco[i] = disco(x, y.d, R = B)$p
  cat('\n  Done. \n')
}
RES = rbind(res.slice, res.disco)
RES.dat = as.data.frame(RES)
colnames(RES.dat) = c(2, 5, 10, 23, 46, 115)
rownames(RES.dat) = c('ECCFIC', 'DISCO')
RES.dat

save.image("D:/Desktop/JASAIndependenceTest/results/Example 4.RData")
library(stargazer)
stargazer(RES.dat, summary = F, digits = 4, table.placement = "H")

