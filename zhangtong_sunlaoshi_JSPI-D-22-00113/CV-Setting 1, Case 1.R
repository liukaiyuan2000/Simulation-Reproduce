Rcpp::sourceCpp("funcs.cpp")
library(ggplot2)
library(tidyverse)
library(latex2exp)
DGP = function(n, cov.u, delta) {  
  const = (c(1:50) - 0.5) * pi  
  g = do.call(cbind, lapply(const^(-2), \(x) rnorm(n, 0, sqrt(x))))  
  phi = sqrt(2) * sin(const %*% t(tm))
  x = g %*% phi
  z1 = rnorm(n, -1, 1)
  z2 = runif(n, 0, pi)
  z = cbind(z1, z2)
  u = matrix(rnorm(n*2, 0, sqrt(cov.u)), n, 2)
  w = z + u
  y = rep(0, n)
  error = rnorm(n, 0, 0.5)
  for (i in c(1:n)) {  
    y[i] = g[i, 1] + 3 * g[i, 2] + z[i, ] %*% beta + error[i] + delta * sqrt(inpro(x[i, ], x[i, ], tm)) * abs(z[i, 1] - z[i, 2]) * 1.5  
  }
  return(list(Y = y, X = x, Z = z, W = w))
}
##CV PROCODURE
tgridnum = 100
tm = c(1:tgridnum)*0.01 - 0.005
beta = c(-1.5,1)
sigma = 0.5
n = 100
cov.u = 0.8
delta = 0
rep = 1000
t1 <- Sys.time()
cvresult = cross_validation(rep, n, cov.u, delta, tm)
t2 <- Sys.time()
print(t2-t1)

order(cvresult[,10])
cvresult[160,]
plot(colMeans(cvresult)/n)

dat <- tibble(
  label = 1:10, 
  meanv = colMeans(cvresult)/n, 
  lower = meanv - apply(cvresult/n, 2, sd), 
  upper = meanv + apply(cvresult/n, 2, sd), 
)
fig <- ggplot(data = dat) + 
  geom_point(aes(x = label, y = meanv), size = 3) + 
  geom_errorbar(
    aes(x = label, ymin = lower, ymax = upper), 
    width = 0.6, linewidth = 1.05
  ) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(
    x = "m", y = TeX('mean value of \\mathbf{CV}(m)'), 
    title = TeX('Setting 1, $\\Sigma_1$, Case 1')
  )
fig



