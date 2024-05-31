rm(list=ls())
# Jiang Du
# May 19 2023
# Beijing University of Technology
#Testing for complete independence in high dimensions
#Author: JAMES R. SCHOTT

Gendata <- function(n,m,rho=0){ 
  if(rho==0){ 
    X <- matrix(rnorm(n*m),n,m) 
  }else{ 
    library(MASS)
    Sigma <- matrix(rho,m,m) 
    diag(Sigma) <- 1 
    X <- mvrnorm(n, rep(0,m),Sigma) 
  }
  X 
}

Wnm <- function(X){ 
  N <- nrow(X) 
  n <- N-1
  m <- ncol(X)  
  R1 <- cor(X)
  if(m<=n){ 
    Tn <- -log(det(R1))*(n-(2*m+5)/6) 
  }else{
    Tn <- 0   
  }
  Tn <- max(Tn,1/m/n/m/n)
  df=m*(m-1)/2
  pvalue <- 1-pchisq(Tn, df=df)
  out=list(Tn=Tn,pvalue=pvalue,df=df)
  out 
}

Schottest <- function(X){
  N <- nrow(X)
  n <- N-1
  m <- ncol(X)
  R1 <- cor(X)
  diag(R1) <- 0 
  Tn <- sum(R1^2)/2-m*(m-1)/2/n
  esd <- sqrt(m*(m-1)*(n-1)/(n+2)/n/n)   
  Tn <- Tn/esd
  out=list(Tn=Tn,pvalue=2*(1-pnorm(abs(Tn))))
  out
}


t1 <- Sys.time()  
ltable=vector("list",4) 
alpha <- 0.05
Sim <- 8000 
for(cs in 1:4){
  if((cs/2)<=1){rho=0}else{ 
    rho <- 0.1 
  }
  Ns <- 2^(2:8)+1 
  Ms <- 2^(2:8)
  RES <- NULL  
  for(n in Ns){
    Res <- NULL
    for(m in Ms){
      res <- rep(0,Sim) 
      for(sim in 1:Sim){ 
        cat(n,m,sim,"\r")
        X <- Gendata(n,m,rho)
        if(cs%%2==1){ 
          fit <- Schottest(X) 
        }else{
          fit <- Wnm(X) 
        }
        res[sim] <- 1*(fit$pvalue<alpha)
      }
      Res <- c(Res,mean(res)) 
    }
    RES <- cbind(RES,Res) 
    
  }
  colnames(RES) <- paste("n=",Ns-1)
  rownames(RES) <- paste("m=",Ms)
  ltable[[cs]]=RES
}



t2 <- Sys.time()

print(ltable)
print(t2-t1)