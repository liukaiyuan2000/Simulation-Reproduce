library(bivpois)
n <- 10000
theta.vec = c(1, 1, 2)
xy <- rbp(n, theta.vec)
x <- xy[, 1]
y <- xy[, 2]
s = pmin(x, y)
Sk = cumsum(x)
sn_star = pmin(s, c(0, Sk[-n]))

vn_vec <- function(r){
  res = 1 / gamma(r + 1) / gamma(x - r + 1) / gamma(y - r + 1)
  return(res)
}

cn_vec <- function(n, k, r){
  lower = pmax(0, k - sn_star)
  upper = pmin(k, sn_star)
  if(n == 1){
    return(vn_vec(r)[1])
  } else{
    res = 0
    for(r in lower[n]:upper[n]){
      res = res + cn_vec(n - 1, k - r, r) * vn_vec(r)[n]
    }
    return(res)
  }
}
cn_vec(1, 0, 0)




