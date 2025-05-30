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

# 计算线性模型统计量的函数
linear_sta <- function(x, y, xnew = lx){
  n = length(y)
  p = ncol(x)
  beta_n <- drop(solve(t(x)%*%x) %*% (t(x)%*%y))
  # n*1 vec
  xnew1 <- drop(x %*% beta_n)
  
  x0 <- quantile(drop(xnew1), 0.99)
  residu <- drop(y - xnew)
  residu1 <- drop(y - xnew1)
  
  
  ss <- outer(xnew, xnew, '<=')
  
  # sigma_n^2 --vector with dimensions n*1(sigma is a constant)
  sign2 <- mean(residu1^2)
  
  # Vn --vector with dimensions n*1
  Vn <- colSums(ss * residu1) / sqrt(n)
  
  # An --vector with dimension n*1
  alpha <- xnew / sign2
  An <- colMeans((alpha^2) * sign2 * t(ss))
  
  # Tnrn1 --vector with dimensions n*1
  sum1 <- rep(0, n)
  for(i in 1:n){
    temp1 <- matrix(ss[, i] * xnew / An, 1, n)
    temp2 <- matrix(residu1 * xnew, n, 1)
    sum1[i] <- temp1 %*% ss %*% temp2 / sign2 / sqrt(n) / n
  }
  Tnrn1 <- Vn - sum1
  
  ## test statistic(Wn) computation
  Wn <- mean((Tnrn1^2) * (xnew1 <= x0)) / 
    mean(xnew1 <= x0) / mean(xnew1 <= x0) / sign2
  return(Wn)
}


# 结果输出
linear_Wn <- linear_sta(x, y)
linear_pvalue <- mean(1 * (linear_Wn <= TN))  # 0
cat('Wn = ', round(linear_Wn, 4), '\n', 'P_value = ', linear_pvalue, sep = "")

# 图
library(ggplot2)
library(latex2exp)
# TeX('$\\hat{\\beta}_n^T X$')
p <- ggplot(data = data.frame(cbind(lx, y)), aes(x = lx, y = y)) + 
  geom_point(color = "firebrick") + xlab(expression(hat(beta)[n]^{T}*X)) + 
  ylab('Y') + xlim(c(-3.3, 2)) + ylim(c(-4, 3)) + 
  theme(panel.grid = element_blank())
p







