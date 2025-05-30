rm(list = ls())
library(kernlab)
library(energy)
library(ggplot2)
library(reshape)
library(ggpubr)
source("D:/Desktop/JASAIndependenceTest/8014586/source.R")
source("D:/Desktop/JASAIndependenceTest/code/DGP.R")
Rcpp::sourceCpp("D:/Desktop/JASAIndependenceTest/code/function.cpp")
B = 50
nsim = 1000
n.choose = seq(20, 50, by = 5)
RES = NULL
for(n in n.choose){
  res.kernel = rep(0, nsim)
  res.dcov   = rep(0, nsim)
  for(i in 1:nsim){
    cat('Example 6: ', 'n = ', n, '; ', i, '\r', sep = "")
    dat = DGP(example = 6, n = n, p = 5)
    x = dat$x
    y = dat$y
    res.kernel[i] = ptest_kem_rcpp(x, y, B = B)$pvalue
    res.dcov[i] = dcov.test(x, y, R = B)$p
  }
  cat('\n  Done. \n')
  temp = cbind(res.kernel, res.dcov)
  RES = rbind(RES, colMeans(temp <= 0.1))
}
colnames(RES) = c('Kernel', 'DCOV')
RES.dat = as.data.frame(RES)
datlong = melt(RES.dat, id.vars = NULL)
fig4 = ggplot(
  datlong, 
  aes(x = rep(n.choose, 2), y = value, group = variable)
) + geom_line(
  aes(colour = variable, linetype = variable, linewidth = variable), 
  lineend = "round"
) + scale_colour_manual(
  values = c(
    Kernel = 'red', DCOV = 'orange'
  )
) + scale_linetype_manual(
  values = c(
    Kernel = 'dashed', DCOV = 'dotted'
  )
) + scale_linewidth_manual(
  values = c(
    Kernel = 1.5, DCOV = 1.5
  )
) + geom_point(
  aes(shape = variable, size = variable, colour = variable)
) + scale_size_manual(
  values = c(
    Kernel = 5, DCOV = 5
  )
) + labs(
  x = "size", y = "power"
) + theme_classic2() + theme(
  legend.title = element_blank(), 
  legend.text = element_text(size = 20), 
  legend.key.size = unit(1.5, 'cm'), 
  legend.key.width = unit(3, 'cm'), 
  legend.position = "inside", 
  legend.position.inside = c(0.4, 0.5), 
  axis.title = element_text(
    family = 'serif', size = 15
  ), 
  axis.text = element_text(size = 12)
)
fig4

save.image("D:/Desktop/JASAIndependenceTest/results/Example 6.RData")

library(showtext)
setEPS()
postscript("D:/Desktop/JASAIndependenceTest/results/Figure4.eps", width = 8, height = 4)
showtext_begin()
fig4
showtext_end()
dev.off()