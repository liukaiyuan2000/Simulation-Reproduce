# A rank-based adaptive independence test for high-dimensional data
setwd("D:/Desktop/RAT")
source("parallel_funcs.R")
set.seed(123456)

t1 <- Sys.time()
N <- c(50, 100, 200)
P <- c(50, 100, 200, 400)
example.choose = 1:3
case.choose = c("a", "b", "c")

idx = 0
RES.sum = list()
RES.name = c()
cl <- makeSOCKcluster(no_cores)
registerDoSNOW(cl)

cv <- NULL
for(n in N){
  for(p in P){
    cat(
      "Calculate critical value: ", "n = ", n, "; p = ", p, "\n", sep = ""
    )
    ECP <- critical.value(n, p, Sim.cv)
    
    temp1 <- c(quantile(ECP[,1], c(alpha/2,1-alpha/2)),quantile(ECP[,2], c(alpha/2,1-alpha/2)),quantile(ECP[,3], c(alpha/2,1-alpha/2)))
    temp2 <- c(quantile(ECP[,4],1-alpha),quantile(ECP[,5],1-alpha),quantile(ECP[,6],1-alpha),quantile(ECP[,7],alpha))
    cv <- rbind(cv, c(temp1,temp2))
    cat("\nDone ~~~\n\n")
  }
}
write.table(cv, file = "cv2.txt")

for(example in example.choose){
  for(case in case.choose){
    idx = idx + 1
    RES <- NULL
    m = 0
    for(n in N){
      REs <- NULL
      for(p in P){
        cat(
          "Example: ", switch(example, "(1)", "(2)", "(3)"), 
          "; Case: ", case, "; n = ", n, "; p = ", p, "\n", sep = ""
        )
        m = m + 1
        cvs <- cv[m, ]
        ## Simulation result
        RES.store <- Sim_func(cvs, n, p, example, case, Sim)
        REs <- rbind(REs,colMeans(RES.store))
        cat("\nDone ~~~\n\n")
      }
      rownames(REs) <- paste("p = ", P, sep = "")
      RES <- rbind(RES,REs)
    }
    colnames(RES) <- c("sumrho","sumtau","sumfoot","maxrho","maxtau","maxfoot","summax")
    RES.sum[[idx]] = RES
    RES.name = c(
      RES.name, 
      paste(
        "Example: ", switch(example, "(1)", "(2)", "(3)"), 
        "; Case: ", case, sep = ""
      )
    )
  }
}

stopCluster(cl)
names(RES.sum) = RES.name

t2 <- Sys.time()
print((t2 - t1))
print(RES.sum)

# write.csv(RES,file = "C:/Users/Administrator/Documents/Xiangyu Shi/CSSC/dnorm1.csv")




