#' @title Data Generate Process
#'
#' @param n sample size
#' @param p dimension
#' @param example example: 1, 2, 3
#' @param case case: "a", "b", "c"
#'
#' @return the data \eqn{X}.
#' @export
#'
#' @examples n = 50; p = 10; example = 1; case = "a"
#' DGP(n, p, example, case)
DGP <- function(n, p, example, case){
  switch(
    example,
    {
      if(case == "a"){
        X <- matrix(rnorm(n*p),n,p)
      } else if(case == "b") {
        X <- matrix(rcauchy(n*p,0,1),n,p)
      } else if(case == "c") {
        X1 <- rnorm(n*p)
        index <- sample(1:(n*p),0.05*n*p)
        X2<- rnorm(0.05*n*p,0,10)
        X1[index] <- X1[index]+X2
        X <- matrix(X1,n,p)
      }
    },
    {
      if(case == "a"){
        rho=0.03
        Sigmma <- (1-rho)*diag(p)+rho*matrix(1,p,p)
        X <- mvrnorm(n, rep(0,p),Sigmma)
      } else if(case == "b") {
        Y<-matrix(0,p,n)
        Z<-matrix(rcauchy(n*p,0,1),n,p)
        for (k in 1:p) {
          Y[k,]<-rowSums(Z[,-k])
        }
        X<-Z+t(Y)/10/p
      } else if(case == "c") {
        rho=0.03
        Sigmma <- (1-rho)*diag(p)+rho*matrix(1,p,p)
        X1 <- mvrnorm(n, rep(0,p),Sigmma)
        X1<-c(X1)
        index <- sample(1:(n*p),0.05*n*p)
        X2<- rnorm(0.05*n*p,0,10)
        X1[index] <- X1[index]+X2
        X <- matrix(X1,n,p)
      }
    },
    {
      if(case == "a"){
        sigmma=matrix(0,p,p)

        Sigmma <- 1*diag(p)
        Sigmma[2,1]<-2.5*sqrt(log(p)/n)
        Sigmma[1,2]<-2.5*sqrt(log(p)/n)

        X <- mvrnorm(n, rep(0,p),Sigmma)
      } else if(case == "b") {
        Z=matrix(rcauchy(n*p,0,1),n,p)
        X=cbind(Z[,1]+sqrt(log(p)/n)*Z[,2],Z[,2]+sqrt(log(p)/n)*Z[,1],Z[,(3:p)])
      } else if(case == "c") {
        sigmma=matrix(0,p,p)

        Sigmma <- 1*diag(p)
        Sigmma[2,1]<-2.5*sqrt(log(p)/n)
        Sigmma[1,2]<-2.5*sqrt(log(p)/n)
        X1 <- mvrnorm(n, rep(0,p),Sigmma)
        X1<-c(X1)
        index <- sample(1:(n*p),0.05*n*p)
        X2<- rnorm(0.05*n*p,0,10)
        X1[index] <- X1[index]+X2
        X <- matrix(X1,n,p)
      }
    }
  )
  return(X)
}


