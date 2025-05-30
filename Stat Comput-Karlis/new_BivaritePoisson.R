m <- 10
theta1 <- 1
theta2 <- 1
X <- rpois(m,theta1)
Y <- rpois(m,theta2)

Sn <- sum(pmin(X,Y))
sn_star <- function(n) min(min(X[n],Y[n]),sum(pmin(X[1:n],Y[1:n])))
vn_vec <- function(n,r){
  #r=先验分布的成分数-1
  vec <- rep(0,Sn+r+1)
  for (i in 0:(Sn+r)) {
    if(X[n]-i>=0&Y[n]-i>=0){
      vec[i+1] <- 1/gamma(X[n]-i+1)/gamma(Y[n]-i+1)/gamma(i+1)
    }else{
      vec[i+1] <- 0
    }
  }
  return(vec)
}
# vn_vec(n=50,r=1)
cn_vec <- function(n,k,r){
  vec <- rep(0,Sn+r+1)
  for (i in 1:n) {
    cat(i, '\r')
    if(i==1){
      vec <- vn_vec(i,r)
    }else{
      lower <- max(0,k-sn_star(i))
      upper <- min(k,sn_star(i))
      for (j in lower:upper) {
        vec[j] <- vn_vec(i,r)[j+1]*cn_vec(i-1,k,r)[k-j+1]
      }
    }
  }
  return(sum(vec))
}
cn_vec(10,Sn,1)


