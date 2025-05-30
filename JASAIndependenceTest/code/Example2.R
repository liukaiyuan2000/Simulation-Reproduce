rm(list = ls())
library(kernlab)
library(energy)
library(ggplot2)
library(ggpubr)
library(reshape2)
source("D:/Desktop/JASAIndependenceTest/8014586/source.R")
source("D:/Desktop/JASAIndependenceTest/code/DGP.R")
Rcpp::sourceCpp("D:/Desktop/JASAIndependenceTest/code/function.cpp")
B = 50
nsim = 1000
p.choose = c(1, seq(5, 100, by = 5))
delta.choose = seq(0, 2, by = 0.2)
RES.power = NULL
for(model.choose in 1:2){
  for(p in p.choose){
    res.slice  = rep(0, nsim)
    res.dcov   = rep(0, nsim)
    res.disco  = rep(0, nsim)
    for(i in 1:nsim){
      cat('Example 2: ', 'p = ', p, '; model = ', model.choose, '; ', i, '\r', sep = "")
      dat = DGP(example = 2, model = model.choose, delta = 1, p = p)
      x = dat$x
      y.d = dat$y.d
      # Discrete y
      res.slice[i] = ptest_slm_rcpp(x, y.d, B)$pvalue
      res.dcov[i] = dcov.test(x, y.d, R = B)$p
      res.disco[i] = disco(x, y.d, R = B)$p
    }
    cat('\n  Done. \n')
    temp = cbind(res.slice, res.dcov, res.disco)
    RES.power = rbind(RES.power, colMeans(temp <= 0.1))
  }
}
RES.delta = NULL
for(model.choose in 1:2){
  for(delta in delta.choose){
    res.slice  = rep(0, nsim)
    res.dcov   = rep(0, nsim)
    res.disco  = rep(0, nsim)
    for(i in 1:nsim){
      cat('Example 2: ', 'delta = ', delta, '; model = ', model.choose, '; ', i, '\r', sep = "")
      dat = DGP(example = 2, model = model.choose, delta = delta, p = 10)
      x = dat$x
      y.d = dat$y.d
      res.slice[i] = ptest_slm_rcpp(x, y.d, B)$pvalue
      res.dcov[i] = dcov.test(x, y.d, R = B)$p
      res.disco[i] = disco(x, y.d, R = B)$p
    }
    cat('\n  Done. \n')
    temp = cbind(res.slice, res.dcov, res.disco)
    RES.delta = rbind(RES.delta, colMeans(temp <= 0.1))
  }
}
### Plots ----
{
data.figb = as.data.frame(RES.power[1:21, ])
data.figa = as.data.frame(RES.delta[1:11, ])
data.figd = as.data.frame(RES.power[22:42, ])
data.figc = as.data.frame(RES.delta[12:22, ])
colnames(data.figa) = c('ECCFIC', 'DCOV', 'DISCO')
colnames(data.figb) = c('ECCFIC', 'DCOV', 'DISCO')
colnames(data.figc) = c('ECCFIC', 'DCOV', 'DISCO')
colnames(data.figd) = c('ECCFIC', 'DCOV', 'DISCO')

datalong.a = melt(data.figa, id.vars = NULL)
datalong.b = melt(data.figb, id.vars = NULL)
datalong.c = melt(data.figc, id.vars = NULL)
datalong.d = melt(data.figd, id.vars = NULL)

figa = ggplot(
  datalong.a, 
  aes(x = rep(delta.choose, 3), y = value, group = variable)
) + geom_line(
  aes(colour = variable, linetype = variable, linewidth = variable), 
  lineend = "round"
) + scale_colour_manual(
  values = c(
    ECCFIC = 'blue', DCOV = 'orange', DISCO = 'magenta'
  )
) + scale_linetype_manual(
  values = c(
    ECCFIC = 'solid', DCOV = 'dotted', DISCO = 'dotdash'
  )
) + scale_linewidth_manual(
  values = c(
    ECCFIC = 1.5, DCOV = 1.5, DISCO = 1.5
  )
) + geom_point(
  aes(shape = variable, size = variable, colour = variable)
) + scale_size_manual(
  values = c(
    ECCFIC = 5, DCOV = 5, DISCO = 5
  )
) + labs(
  x = expression(delta), y = "power", title = "(a)"
) + theme_classic2() + theme(
  legend.title = element_blank(), 
  legend.text = element_text(size = 20), 
  legend.key.size = unit(1.5, 'cm'), 
  legend.key.width = unit(3, 'cm'), 
  legend.position = "none", 
  legend.position.inside = c(0.8, 0.4), 
  axis.title = element_text(
    family = 'serif', size = 18
  ), 
  axis.text = element_text(size = 15), 
  plot.title = element_text(colour = 'red')
)
figa

figb = ggplot(
  datalong.b, 
  aes(x = rep(p.choose, 3), y = value, group = variable)
) + geom_line(
  aes(colour = variable, linetype = variable, linewidth = variable), 
  lineend = "round"
) + scale_colour_manual(
  values = c(
    ECCFIC = 'blue', DCOV = 'orange', DISCO = 'magenta'
  )
) + scale_linetype_manual(
  values = c(
    ECCFIC = 'solid', DCOV = 'dotted', DISCO = 'dotdash'
  )
) + scale_linewidth_manual(
  values = c(
    ECCFIC = 1.5, DCOV = 1.5, DISCO = 1.5
  )
) + geom_point(
  aes(shape = variable, size = variable, colour = variable)
) + scale_size_manual(
  values = c(
    ECCFIC = 5, DCOV = 5, DISCO = 5
  )
) + labs(
  x = 'p', y = "power", title = "(b)"
) + theme_classic2() + theme(
  legend.title = element_blank(), 
  legend.text = element_text(size = 20), 
  legend.key.size = unit(1.5, 'cm'), 
  legend.key.width = unit(3, 'cm'), 
  legend.position = "none", 
  legend.position.inside = c(0.8, 0.4), 
  axis.title = element_text(
    family = 'serif', size = 18
  ), 
  axis.text = element_text(size = 15), 
  plot.title = element_text(colour = 'red')
)
figb

figc = ggplot(
  datalong.c, 
  aes(x = rep(delta.choose, 3), y = value, group = variable)
) + geom_line(
  aes(colour = variable, linetype = variable, linewidth = variable), 
  lineend = "round"
) + scale_colour_manual(
  values = c(
    ECCFIC = 'blue', DCOV = 'orange', DISCO = 'magenta'
  )
) + scale_linetype_manual(
  values = c(
    ECCFIC = 'solid', DCOV = 'dotted', DISCO = 'dotdash'
  )
) + scale_linewidth_manual(
  values = c(
    ECCFIC = 1.5, DCOV = 1.5, DISCO = 1.5
  )
) + geom_point(
  aes(shape = variable, size = variable, colour = variable)
) + scale_size_manual(
  values = c(
    ECCFIC = 5, DCOV = 5, DISCO = 5
  )
) + labs(
  x = expression(delta), y = "power", title = "(c)"
) + theme_classic2() + theme(
  legend.title = element_blank(), 
  legend.text = element_text(size = 20), 
  legend.key.size = unit(1.5, 'cm'), 
  legend.key.width = unit(3, 'cm'), 
  legend.position = "none", 
  axis.title = element_text(
    family = 'serif', size = 18
  ), 
  axis.text = element_text(size = 15), 
  plot.title = element_text(colour = 'red')
)
figc

figd = ggplot(
  datalong.d, 
  aes(x = rep(p.choose, 3), y = value, group = variable)
) + geom_line(
  aes(colour = variable, linetype = variable, linewidth = variable), 
  lineend = "round"
) + scale_colour_manual(
  values = c(
    ECCFIC = 'blue', DCOV = 'orange', DISCO = 'magenta'
  )
) + scale_linetype_manual(
  values = c(
    ECCFIC = 'solid', DCOV = 'dotted', DISCO = 'dotdash'
  )
) + scale_linewidth_manual(
  values = c(
    ECCFIC = 1.5, DCOV = 1.5, DISCO = 1.5
  )
) + geom_point(
  aes(shape = variable, size = variable, colour = variable)
) + scale_size_manual(
  values = c(
    ECCFIC = 5, DCOV = 5, DISCO = 5
  )
) + labs(
  x = 'p', y = "power", title = "(d)"
) + theme_classic2() + theme(
  legend.title = element_blank(), 
  legend.text = element_text(size = 20), 
  legend.key.size = unit(1.5, 'cm'), 
  legend.key.width = unit(3, 'cm'), 
  legend.position = "none", 
  axis.title = element_text(
    family = 'serif', size = 18
  ), 
  axis.text = element_text(size = 15), 
  plot.title = element_text(colour = 'red')
)
figd
fig1 = ggarrange(
  figa, figb, figc, figd, nrow = 2, ncol = 2, common.legend = T
)
fig1
}
#####

save.image("D:/Desktop/JASAIndependenceTest/results/Example 2.RData")

library(showtext)
setEPS()
postscript("D:/Desktop/JASAIndependenceTest/results/Figure1.eps", width = 12, height = 10)
showtext_begin()
fig1
showtext_end()
dev.off()












