rm(list = ls())
library(kernlab)
library(energy)
library(ggplot2)
library(reshape)
library(ggpubr)
source("D:/Desktop/JASAIndependenceTest/8014586/source.R")
source("D:/Desktop/JASAIndependenceTest/code/DGP.R")
B = 50
nsim = 1000
RES = NULL
for(n in seq(10, 50, by = 5)){
  res.slice  = rep(0, nsim)
  res.kernel = rep(0, nsim)
  res.dcov   = rep(0, nsim)
  res.disco  = rep(0, nsim)
  for(i in 1:nsim){
    cat('Example 5: ', 'n = ', n, '; ', i, '\r', sep = "")
    dat = DGP(example = 5, n = n)
    x = dat$x
    y = dat$y
    y.d = dat$y.d
    # Discrete y
    ptest.s.res <- ptest_H2d(x, y.d, B = B)
    ptest.k.res <- ptest_H2c(x, y, B = B)
    res.slice[i] = ptest.s.res$pvalue
    res.kernel[i] = ptest.k.res$pvalue
    res.dcov[i] = dcov.test(x, y, R = B)$p
    res.disco[i] = disco(x, y.d, R = B)$p
  }
  cat('\n  Done. \n')
  temp = cbind(res.slice, res.kernel, res.dcov, res.disco)
  RES = rbind(RES, colMeans(temp <= 0.1))
}
colnames(RES) = c('Slice', 'Kernel', 'DCOV', 'DISCO')
RES.dat = as.data.frame(RES)
datlong = melt(RES.dat, id.vars = NULL)
fig3 = ggplot(
  datlong, 
  aes(x = rep(seq(10, 50, by = 5), 4), y = value, group = variable)
) + geom_line(
  aes(colour = variable, linetype = variable, linewidth = variable), 
  lineend = "round"
) + scale_colour_manual(
  values = c(
    Slice = 'blue', Kernel = 'red', 
    DCOV = 'orange', DISCO = 'magenta'
  )
) + scale_linetype_manual(
  values = c(
    Slice = 'solid', Kernel = 'dashed', 
    DCOV = 'dotted', DISCO = 'dotdash'
  )
) + scale_linewidth_manual(
  values = c(
    Slice = 1.5, Kernel = 1.5, 
    DCOV = 1.5, DISCO = 1.5
  )
) + geom_point(
  aes(shape = variable, size = variable, colour = variable)
) + scale_size_manual(
  values = c(
    Slice = 5, Kernel = 5, 
    DCOV = 5, DISCO = 5
  )
) + labs(
  x = "size", y = "power"
) + theme_classic2() + theme(
  legend.title = element_blank(), 
  legend.text = element_text(size = 20), 
  legend.key.size = unit(1.5, 'cm'), 
  legend.key.width = unit(3, 'cm'), 
  legend.position = "inside", 
  legend.position.inside = c(0.8, 0.4), 
  axis.title = element_text(
    family = 'serif', size = 15
  ), 
  axis.text = element_text(size = 12)
)
fig3

save.image("D:/Desktop/JASAIndependenceTest/results/Example 5.RData")

library(showtext)
setEPS()
postscript("D:/Desktop/JASAIndependenceTest/results/Figure3.eps", width = 8, height = 4)
showtext_begin()
fig3
showtext_end()
dev.off()