rm(list = ls())
Rcpp::sourceCpp("D:/Desktop/非参数核带宽选择/cpp_funcs.cpp")
## Bandwidth Selection in Nonparametric Kernel Testing - JASA
## Only Proposed Method
n = 50
B = 250
nsim = 500
K_square <- function(x){
  dnorm(x)^2
}
int_K_square <- integrate(K_square, lower = -Inf, upper = Inf)$value
## calculate K^{(3)}(0) in \hat{a}_1
Kern_conv3_zero <- 1/(4*sqrt(2*pi))

Sim_cpp(
  n, B, nsim, 
  Model_case = 1, H_case = 2, j_case = 1, k = 1, 
  alpha = 0.99, int_K_square, Kern_conv3_zero, 
  set_h = seq(0.01, 0.99, 0.01)
)




