rm(list = ls())
source("D:/桌面/ICM模拟/Tnsztest/Tnsz_BMfn().R")

# 数据读取
mydata <- read.table('D:/桌面/杜老师论文220918/airfoil_self_noise.dat')
names(mydata) <- c('X1', 'X2', 'X3', 'X4', 'X5', 'Y')
mydata <- as.matrix(mydata)

# 标准化
mydata2 <- scale(mydata)
beta1 <- c(-0.6323, -0.4339, -0.5339, 0.2386, -0.2644)
x <- mydata2[, -6]
y <- mydata2[, 6]
lx <- drop(x %*% beta1)

# 生成一组标准Brown运动随机数方便计算P_value
Sim <-  5000
TN <- rep(0, Sim)
for(sim in 1:Sim){
  cat(sim,'\r')
  X <- BMfn(M = 500, T1 = 1)
  TN[sim] <- mean(X^2)
}

model <- lm(y ~ lx + I(lx^2) + I(lx^3))
theta_n <- model$coef

polynomial_sta <- function(x, y, theta = theta_n){
  # model: m(x) = x + x^2 + x^3
  n = length(y)
  p = ncol(x)
  beta_n <- drop(solve(t(x)%*%x) %*% (t(x)%*%y))
  beta_n <- beta_n / sqrt(sum(beta_n^2))
  # beta_n <- beta1
  B_n <- beta1
  alpha_hat <- 1
  # n*1 vec
  xnew <- drop(x %*% B_n %*% alpha_hat)
  xnew1 <- drop(x %*% beta_n)
  
  x0 <- quantile(drop(xnew1), 0.99)
  residu1 <- drop(y - theta[1] - theta[2] * xnew1 - 
                    theta[3] * xnew1^2 - theta[4] * xnew1^3)
  
  
  ss <- outer(xnew, xnew, '<=')
  
  # sigma_n^2 --vector with dimensions n*1(sigma is a constant)
  sign2 <- mean(residu1^2)
  
  # Vn --vector with dimensions n*1
  Vn <- colSums(ss * residu1) / sqrt(n)
  
  # An --vector with dimension n*1
  dm <- theta_n[2] + 2 * theta_n[3] * xnew + 3 * theta_n[4] * xnew^2
  alpha <- dm * xnew / sign2
  An <- colMeans((alpha^2) * sign2 * t(ss))
  
  # Tnrn1 --vector with dimensions n*1
  sum1 <- rep(0, n)
  for(i in 1:n){
    temp1 <- matrix(ss[, i] * xnew * dm / An, 1, n)
    temp2 <- matrix(residu1 * alpha, n, 1)
    sum1[i] <- temp1 %*% ss %*% temp2 / sqrt(n) / n
  }
  Tnrn1 <- Vn - sum1
  
  ## test statistic(Wn) computation
  Wn <- mean((Tnrn1^2) * (xnew1 <= x0)) / 
    mean(xnew1 <= x0) / mean(xnew1 <= x0) / sign2
  return(Wn)
}


# 结果输出
polynomial_Wn <- polynomial_sta(x, y)
polynomial_pvalue <- mean(1 * (polynomial_Wn <= TN))  # 0.2414
cat('Wn = ', round(polynomial_Wn, 4), '\n', 
    'P_value = ', polynomial_pvalue, '\n', sep = "")


# 图
library(ggplot2)
library(latex2exp)
# TeX('$\\hat{\\beta}_n^T X$')
p <- ggplot(data = data.frame(cbind(lx, y)), aes(x = lx, y = y)) + 
  geom_point(color = "firebrick") + xlab(expression(hat(beta)[n]^{T}*X)) + 
  ylab('Y') + xlim(c(-3.3, 2)) + ylim(c(-4, 3)) + 
  theme(panel.grid = element_blank()) + 
  geom_smooth(method="lm", formula = y ~ x + I(x^2) + I(x^3), se=FALSE)
p








