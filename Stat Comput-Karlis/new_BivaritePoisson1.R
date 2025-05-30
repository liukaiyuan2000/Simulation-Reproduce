n <- 50
theta1 <- 1
theta2 <- 1
x <- rpois(n,theta1)
y <- rpois(n,theta2)
s = pmin(x, y)
Sk = cumsum(x)
sn_star = pmin(s, c(0, Sk[-n]))

vn_vec <- function(r){
  res = 1 / gamma(r + 1) / gamma(x - r + 1) / gamma(y - r + 1)
  return(res)
}

cn_vec <- function(k, r){
  lower = pmax(0, k - sn_star)
  upper = pmin(k, sn_star)
  res = rep(0, n)
  res[1] = vn_vec(k)[1]
  for(i in 2:n){
    for(r in lower:upper){
      res[i] = res[i] + cn_vec(k - r, r)[i - 1] * vn_vec(r)[i]
    }
  }
  
  return(res)
}
cn_vec(50,Sn,9)



