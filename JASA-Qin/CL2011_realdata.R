source("D:/Desktop/JASA-Qin/CL2011_func.R")
library(plsgenomics)
library(corrplot)
data(SRBCT)
X = SRBCT$X[1:63, ]
Y = SRBCT$Y[1:63]
F_values = cal_F_statistics(X, Y)
top_40_gene_idx = order(F_values, decreasing = T)[1:40]
bottom_160_gene_idx = order(F_values, decreasing = F)[1:160]

X.new = cbind(X[, top_40_gene_idx], X[, bottom_160_gene_idx])
Sigma.hat = cov(X.new)
thre.cv.AL = thre_cv_func(X.new, s_type = 'adaptive lasso', delta.len = 100)
thre.cv.Hard = thre_cv_func(X.new, s_type = 'hard', delta.len = 100)
thre.cv.AL.u = thre_cv_func(X.new, s_type = 'adaptive lasso', delta.len = 100, method = 'universal')
thre.cv.Hard.u = thre_cv_func(X.new, s_type = 'hard', delta.len = 100, method = 'universal')
Sigma.AL = thresholding(Sigma.hat, thre.cv.AL)
Sigma.Hard = thresholding(Sigma.hat, thre.cv.Hard)
Sigma.AL.u = thresholding(Sigma.hat, thre.cv.AL.u)
Sigma.Hard.u = thresholding(Sigma.hat, thre.cv.Hard.u)
par(mfrow = c(2, 2))
fig.hard <- corrplot(
  Sigma.Hard > 0, bg = "white", cl.pos = "n", 
  tl.pos = "n", method = "shade", 
  title = paste(
    'Hard: ', 
    round(sum(Sigma.Hard == 0)/200^2*100, 2), 
    '% zeros', sep = ""), 
  mar = c(0,0,2,0)
)
fig.al <- corrplot(
  Sigma.AL > 0, bg = "white", cl.pos = "n", 
  tl.pos = "n", method = "shade", 
  title = paste(
    'Adaptive Lasso: ', 
    round(sum(Sigma.AL == 0)/200^2*100, 2), 
    '% zeros', sep = ""), 
  mar = c(0,0,2,0)
)
fig.hard.u <- corrplot(
  Sigma.Hard.u > 0, bg = "white", cl.pos = "n", 
  tl.pos = "n", method = "shade", 
  title = paste(
    'Hard (Universal): ', 
    round(sum(Sigma.Hard.u == 0)/200^2*100, 2), 
    '% zeros', sep = ""), 
  mar = c(0,0,2,0)
)
fig.al.u <- corrplot(
  Sigma.AL.u > 0, bg = "white", cl.pos = "n", 
  tl.pos = "n", method = "shade", 
  title = paste(
    'Adaptive Lasso (Universal): ', 
    round(sum(Sigma.AL.u == 0)/200^2*100, 2), 
    '% zeros', sep = ""), 
  mar = c(0,0,2,0)
)
par(mfrow = c(1, 1))







