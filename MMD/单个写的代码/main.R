#A Kernel Method for the Two-Sample-Problem 
rm(list = ls())
set.seed(123)
library(MASS)
library(Rcpp)
library(eummd)
t1 <- Sys.time()
n <- 40
m <- 40
d.choose <- c(20, 40, 80, 100)
K.choose <- c(10, 50)

B = 100
Sim = 200
alpha = 0.05
# Need to change in different case
delta = 0; rho = 0
RES = c()
for(d in d.choose){
  for(K in K.choose){
    for (sim in 1:Sim) {
      cat(
        "d = ", d, "; K = ", K, 
        "; (", sim, "/", Sim, ")Replication\r", 
        sep = ""
      )
      X <- matrix(rnorm(n*d, 0, 1), n, d)
      Y <- matrix(rnorm(m*d, 0, 1), m, d)
      
      Tn <- meammd(
        X = X, Y = Y, pval = TRUE, type = "proj", 
        numproj = K, numperm = 0
      )$stat
      ECp <- matrix(0, Sim, B)
      for (b in 1:B) {
        Z <- rbind(X, Y)
        set <- sample(1:(n+m), n)
        a <- as.vector(set)
        
        X1 <- as.matrix(Z[a,])
        Y1 <- as.matrix(Z[-a,])
        
        Tn_b <- meammd(
          X = X1, Y = Y1, pval = TRUE, type = "proj", 
          numproj = K, numperm = 0
        )$stat
        ECp[sim,b] <- 1*(Tn_b >= Tn)
      }
    }
    RES <- cbind(RES, mean(rowMeans(ECp) > alpha))
  }
}

t2 <- Sys.time()
print((t2-t1))
print(RES)




