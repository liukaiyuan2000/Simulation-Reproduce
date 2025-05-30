DGP <- function(n, p, example.choose, case.choose){
  switch(
    example.choose, 
    {
      if(case.choose == "a"){
        X <- matrix(rnorm(n*p),n,p)
      } else if(case.choose == "b") {
        X <- matrix(rcauchy(n*p,0,1),n,p)
      } else if(case.choose == "c") {
        X1 <- rnorm(n*p)
        index <- sample(1:(n*p),0.05*n*p)
        X2<- rnorm(0.05*n*p,0,10)
        X1[index] <- X1[index]+X2
        X <- matrix(X1,n,p)
      }
    }, 
    {
      if(case.choose == "a"){
        rho=0.03
        Sigmma <- (1-rho)*diag(p)+rho*matrix(1,p,p)
        X <- mvrnorm(n, rep(0,p),Sigmma)
      } else if(case.choose == "b") {
        Y<-matrix(0,p,n)
        Z<-matrix(rcauchy(n*p,0,1),n,p)
        for (k in 1:p) {
          Y[k,]<-rowSums(Z[,-k])
        }
        X<-Z+t(Y)/10/p
      } else if(case.choose == "c") {
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
      if(case.choose == "a"){
        sigmma=matrix(0,p,p)

        Sigmma <- 1*diag(p)
        Sigmma[2,1]<-2.5*sqrt(log(p)/n)
        Sigmma[1,2]<-2.5*sqrt(log(p)/n)

        X <- mvrnorm(n, rep(0,p),Sigmma)
      } else if(case.choose == "b") {
        Z=matrix(rcauchy(n*p,0,1),n,p)
        X=cbind(Z[,1]+sqrt(log(p)/n)*Z[,2],Z[,2]+sqrt(log(p)/n)*Z[,1],Z[,(3:p)])
      } else if(case.choose == "c") {
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


