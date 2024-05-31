##### New test statistics and their wild bootstrap counterparts ----
rm(list = ls())
M.tests = function(X,Y,n,epsilon = 0.0001){
  
  F = ecdf(X)
  G = ecdf(Y)
  Z = c(X,Y)
  E = F(Z)-G(Z)
  KS = sqrt(n)*max(abs(E))
  
  H = ecdf(pmax(X,Y))
  w = sqrt( abs( F(Z)+G(Z)-2*H(Z)-(F(Z)-G(Z))^2 ) )
  wKS = sqrt(n)*max(abs(ifelse(w == 0,0,E/w)))
  
  Z.t = sort(Z)[ceiling(epsilon*2*n):floor((1-epsilon)*2*n)]
  E.t = F(Z.t)-G(Z.t)
  w.t = sqrt( abs(F(Z.t)+G(Z.t)-2*H(Z.t)-(F(Z.t)-G(Z.t))^2 ) )
  WKS = sqrt(n)*max(abs(ifelse(w.t == 0,0,E.t/w.t)))
  
  M = max(KS,WKS)
  
  all = c(KS,wKS,WKS,M)
  return(all)
  
}

M.tests.boot = function(X,Y,n,epsilon = 0.0001,vec.Xi, par.c){
  
  Z = c(X,Y)
  Z.s = sort(Z)
  E.boot = matrix(0,1,2*n)
  for(i in 1:(2*n)){
    E.boot[1,i] = sum( vec.Xi*( ifelse(X <= Z.s[i],1,0) - ifelse(Y <= Z.s[i],1,0) ) )/sqrt(n)
  }
  KS.boot = max(abs(E.boot))
  
  w.boot = matrix(0,1,2*n)
  for(i in 1:(2*n)){
    w.boot[1,i] = sum( abs(vec.Xi/sqrt(2/pi))*( ifelse(X <= Z.s[i],1,0) + ifelse(Y <= Z.s[i],1,0) - 2*ifelse(pmax(X,Y) <= Z.s[i],1,0) ) )/n -
      (sum( abs(vec.Xi/sqrt(2/pi))*(ifelse(X <= Z.s[i],1,0) - ifelse(Y <= Z.s[i],1,0)) ))^2/n^2 
  }
  wKS.boot = max(abs(ifelse(w.boot == 0,0,E.boot/sqrt(abs(par.c*w.boot)))))
  
  Z.t = sort(Z)[ceiling(epsilon*2*n):floor((1-epsilon)*2*n)]
  dl = length(Z.t)
  E.t.boot = matrix(0,1,dl)
  w.t.boot = matrix(0,1,dl)
  for(i in 1:dl){
    E.t.boot[1,i] = sum( vec.Xi*( ifelse(X <= Z.t[i],1,0) - ifelse(Y <= Z.t[i],1,0) ) )/sqrt(n)
    w.t.boot[1,i] = sum( abs(vec.Xi/sqrt(2/pi))*( ifelse(X <= Z.t[i],1,0) + ifelse(Y <= Z.t[i],1,0) - 2*ifelse(pmax(X,Y) <= Z.t[i],1,0) ) )/n -
      (sum( abs(vec.Xi/sqrt(2/pi))*(ifelse(X <= Z.t[i],1,0) - ifelse(Y <= Z.t[i],1,0)) ))^2/n^2 
  }
  
  WKS.boot = max(abs(ifelse(w.t.boot == 0,0,E.t.boot/sqrt(abs(par.c*w.t.boot)))))
  
  M.boot = max(KS.boot,WKS.boot)
  
  all.boot = c(KS.boot,wKS.boot,WKS.boot,M.boot)
  return(all.boot)
  
}

###################################################

## initial set ----
# install.packages('MyMtest_0.1.0.tar.gz', repos = NULL, type = 'source')
library(doSNOW)
library(tcltk)
library(MyMtest)
options(warn = -1)
alpha = 0.05

n = 50

time = c(rep(0,n), rep(1,n))  # auxiliary objects
subject = c(1:n,1:n)          # for BMP.test 

hh = c(1/2,1,2,3,4)	# bandwidth window parameters for Martinez-Camblor (2010) test statistic
B0 = 499			# and the related
B1 = 100			# numbers of repetitions

cv.VGH = 3.161	# critical value of the Vexler et al. (2013) test
# under alpha = 0.05 and n = 50, cf. Table 1, p.339, 

par.c = 0.8		# parameter c in the paper

Nk = 1000		# number of Monte Carlo runs 
Nk.boot = 1000	# number of wild bootstrap runs

no_cores <- 8
pb <- txtProgressBar(max = Nk, style = 3, char = "*")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
#####

## Origin code ----
P_KS  = matrix(0, 1, Nk) 
P_M   = matrix(0, 1, Nk) 
valStat = matrix(0, 2, Nk)
start.time <- Sys.time()
for(i in 1:Nk){
  ###################################
  ### FGM copula, Generator begins 
  ###################################
  
  theta = 0.5		### this parameter changes and it occurs in the tables as $theta$
  
  X = runif(n)
  T = runif(n)
  
  A = 1 + theta*(1-2*X)
  B = sqrt( A^2 - 4*(A-1)*T )
  Y = 2*T/(A+B)
  
  ##################################
  ### FGM copula, Generator ends 
  ##################################
  
  Z = c(X,Y)
  
  score = M.tests(X,Y,n,epsilon = 0.0001)
  
  valStat[1,i] = score[1]
  valStat[2,i] = score[4]
  val_M_boot = matrix(0, 2, Nk.boot)
  
  for(j in 1:Nk.boot){
    cat('i = ', i, '; j = ', j, '\r', sep = "")
    Xi = rnorm(n)
    outcome = M.tests.boot(X,Y,n,epsilon = 0.0001,Xi,par.c)
    val_M_boot[1, j] = outcome[1]
    val_M_boot[2, j] = outcome[4]
  }
  cat('\n')
  p.v.KS = mean( valStat[1,i] < val_M_boot[1,] )
  P_KS[1,i] = ifelse( p.v.KS < alpha, 1,0)
  
  p.v.M = mean( valStat[2,i] < val_M_boot[2,] )
  P_M[1,i] = ifelse( p.v.M < alpha, 1,0)
}
end.time <- Sys.time()
difftime(end.time, start.time, units = 'secs')
#####

## Rcpp for bootstrap ----
P_KS  = matrix(0, 1, Nk) 
P_M   = matrix(0, 1, Nk) 
valStat = matrix(0, 2, Nk)
start.time <- Sys.time()
for(i in 1:Nk){
  cat('i = ', i, '\r', sep = "")
  ###################################
  ### FGM copula, Generator begins 
  ###################################
  
  theta = 0.5		### this parameter changes and it occurs in the tables as $theta$
  
  X = runif(n)
  T = runif(n)
  
  A = 1 + theta*(1-2*X)
  B = sqrt( A^2 - 4*(A-1)*T )
  Y = 2*T/(A+B)
  
  ##################################
  ### FGM copula, Generator ends 
  ##################################
  
  Z = c(X,Y)
  
  score = M_tests(X,Y,n,epsilon = 0.0001)
  
  valStat[1,i] = score[1]
  valStat[2,i] = score[4]
  val_M_boot <- M_boot_func(Nk.boot, X, Y, n, epsilon = 0.0001, par.c)
  
  p.v.KS = mean( valStat[1,i] < val_M_boot[1,] )
  P_KS[1,i] = ifelse( p.v.KS < alpha, 1,0)
  
  p.v.M = mean( valStat[2,i] < val_M_boot[2,] )
  P_M[1,i] = ifelse( p.v.M < alpha, 1,0)
}
end.time <- Sys.time()
difftime(end.time, start.time, units = 'secs')
#####

## Rcpp & Parallel ----
cl <- makeSOCKcluster(no_cores)
registerDoSNOW(cl)
start.time <- Sys.time()
valStat <- foreach(
  sim = 1:Nk, .combine = "rbind", .options.snow = opts, 
  .packages = "MyMtest"
) %dopar% {
  ###################################
  ### FGM copula, Generator begins 
  ###################################
  
  theta = 0.5		### this parameter changes and it occurs in the tables as $theta$
  
  X = runif(n)
  T = runif(n)
  
  A = 1 + theta*(1-2*X)
  B = sqrt( A^2 - 4*(A-1)*T )
  Y = 2*T/(A+B)
  
  ##################################
  ### FGM copula, Generator ends 
  ##################################
  
  Z = c(X,Y)
  
  score = M_tests(X,Y,n,epsilon = 0.0001)
  
  val_M_boot <- M_boot_func(Nk.boot, X, Y, n, epsilon = 0.0001, par.c)
  
  p.v.KS = mean( score[1] < val_M_boot[1,] )
  P_KS = ifelse( p.v.KS < alpha, 1,0)
  
  p.v.M = mean( score[4] < val_M_boot[2,] )
  P_M = ifelse( p.v.M < alpha, 1,0)
  return(c(P_KS, P_M))
}
end.time <- Sys.time()
difftime(end.time, start.time, units = 'secs')
close(pb)
stopCluster(cl)
#####

colMeans(valStat)*100
