BMfn <- function(M = 500, T1 = 1){
  #Function generates a standard brownian motion on the [0, T] interval#
  dt <- T1 / (M - 1) #To ensure we start at zero and end at 1
  W <- numeric(M)
  W[1] = 0
  for(j in 2:M){
    W[j] <- W[j-1] + sqrt(dt) * rnorm(1)
  }
  return(W)
}
# print(quantile(TN, c(0.9, 0.95, 0.99)))
