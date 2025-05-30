## ICMsta: input(x, y) x，y -协变量和响应变量
##         output(icm) -统计量
# library(MASS)
ICMsta <- function(x, y){
    beta_n <- drop(solve(t(x)%*%x) %*% (t(x)%*%y))
    residu <- drop(y - x %*% beta_n)
    
    x_new <- drop(x %*% beta_n)
    id1 <- rep(1:n, n)
    id2 <- rep(1:n, each = n)
    X <- x[id1, ] - x[id2, ]
    A <- matrix(apply(X, 1, crossprod), n, n)
    AA <- exp(-0.5 * A)
    resi <- residu %*% t(residu)
    ## proposed test statistic
    icm <- sum(AA * resi) / n
    # icm <- drop(t(residu) %*% AA %*% residu / n)
    # 上一个公式的矩阵运算较慢，需要更改R中的dll文件
    return(icm)
}















