
print(date())

set.seed(1)

##### Wilcoxon test statistic

W.test = function(X,Y,n){

Z = Y - X
Z.a = abs(Z)
R.a = rank(Z.a)
Z.s = sign(Z)
V = n*(n+1)*(2*n+1)/6
WSR = sum(Z.s*R.a)/sqrt(V)
return(WSR)

}

##### Brunner et al. (1999) test statistic is calculated using the package nparLD
##### Noguchi, K., Gel, Y. R., Brunner, E., Konietschke, F. (2012).
##### nparLD: An R Software package for the nonparametric analysis of longitudinal data in factorial experiments.
##### Journal of Statistical Software, 50, 12.

library(MASS)
library(nparLD)

##### Martinez-Camblor (2010) test statistic 

fAC<- function(X,hh,B){k<- length(X[1,]); n<- length(X[,1]);T<- length(hh)
        l<- length(t);ac<- matrix(0,ncol=T,nrow=(B+1));
        XB<- X;
         
        densidades<- function(X,h,k,t)
                      {l<- length(t);D=matrix(0,nrow=l,ncol=k);
                       pr<- t[1];fn<- t[l];
                       for (i in 1:k) {A<-density(X[,i],adjust=h,from=pr,to=fn);D[,i]<-A$y}
                       return(D)}
           
        g<- function(D,t)
             {l<- length(t);ni<-rep(0,l);
   	      for (i in 1:l) {ni[i]=min(D[i,])};
              return((t[2]-t[1])*(sum(ni)-0.5*(ni[1]+ni[l])))};
               
                   
        for (tt in 1:T) {t<- density(c(X),adjust=hh[tt])$x; ac[1,tt]<- g(densidades(X,hh[tt],k,t),t)}; 
             
        for (b in 1:B) 
            {for (j in 1:n) {XB[j,]<- sample(X[j,])};
             BOOT<- NULL; 
             for (i in 1:k) {BOOT<- cbind(BOOT,n^(-1/5)*rnorm(n))}; 
             XB<- XB+BOOT;
             for (tt in 1:T)
                 {t<- density(c(XB),adjust=hh[tt])$x;              
                  DB<- densidades(XB,hh[tt],k,t);             
                  ac[(b+1),tt]<- g(DB,t)}}
       return(ac)}

MC.test <- function(X,hh,B0,B1){

	k<- length(X[1,]); n<- length(X[,1]);T<- length(hh);hbt<- rep(0,B1)
        Pac<- rep(0,T); Pax<- rep(0,T);ac1<- matrix(0,ncol=T,nrow=B1)
                 
        densidades<- function(X,h,k,t)
                      {l<- length(t);D=matrix(0,nrow=k,ncol=l);
                       for (i in 1:k) {hi<- h*var(X[,i])^0.5*n^-0.2;
                                       for (j in 1:l) {D[i,j]<-mean(dnorm(t[j],X[,i],hi))}}
                       return(D)}
           
        g<- function(D,t)
             {l<- length(t);ni<-rep(0,l);
   	      for (i in 1:l) {ni[i]=min(D[,i])};
              return((t[2]-t[1])*(sum(ni)-0.5*(ni[1]+ni[l])))};

        sal<- fAC(X,hh,B0); for (tt in 1:T) {Pac[tt]<- (rank(sal[,tt])[1]-1)/B0}; 
             
        if ( (sum(Pac < 0.05)==0) | (sum(Pac > 0.05)==0) ) {Pj<- min(Pac); Pb<- Pj; hj<- hh[which.min(Pac)]; hb<- hj}else{
              for (b in 1:B1) 
               {idx<- sample(1:n,replace=TRUE);BOOT<- NULL; 
                for (i in 1:k) {BOOT<- cbind(BOOT,n^(-1/5)*rnorm(n))}; 
                XB<- X[idx,]+BOOT;
                for (tt in 1:T) 
                    {t<- density(c(XB),adjust=hh[tt])$x;
                     DB<- densidades(XB,hh[tt],k,t);ac1[b,tt]<- g(DB,t)}}

                for (b in 1:B1) {for (tt in 1:T) {Pax[tt]<- (rank(c(ac1[b,tt],sal[2:(B0+1),tt]))[1]-1)/B0}; hbt[b]<- hh[which.min(Pax)]};
                hb<- mean(hbt); Pb<- (rank(fAC(X,hb,B0))[1]-1)/B0;
                      }
        return(Pb)}

##### Ghattas et al. (2011) test statistic based on Hermite polynomials

 library(copula)
 library(Matrix)
 library(matrixcalc)

hermite_rec_formula = function(j = 1, x){
result = list()
result[["0"]] = 1
result[["1"]] = 2*x
for (i in 1:j){
result[[as.character(i + 1)]] = 2*x*result[[as.character(i)]] -
2*i*result[[as.character(i - 1)]]
}
return(result)
}

GPRY.test = function(sample){

rec_formula = hermite_rec_formula
density_fun_x = dnorm
density_fun_y = dnorm

statistics <- c()
max_j = 15
for (j in 1:max_j){
j_vec = 1:j
Q_tilde_X = sapply(j_vec,function(j){rec_formula(j,sample$X)[[as.character(j)]]}) *density_fun_x(sample$X)
Q_tilde_Y = sapply(j_vec,function(j){rec_formula(j,sample$Y)[[as.character(j)]]}) *density_fun_y(sample$Y)
V_s_k = Q_tilde_X - Q_tilde_Y
U_n_k = apply(V_s_k, 2, sum)/sqrt(nrow(sample))
W_n_k = t(V_s_k)%*%V_s_k/nrow(sample)
matrix_rank = rankMatrix(W_n_k)[1]
if (j > matrix_rank){
j_vec = c(1:matrix_rank, (matrix_rank + 2))
Q_tilde_X = sapply(j_vec,function(j){rec_formula(j,sample$X)[[as.character(j)]]}) *density_fun_x(sample$X)
Q_tilde_Y = sapply(j_vec,function(j){rec_formula(j,sample$Y)[[as.character(j)]]}) *density_fun_y(sample$Y)
V_s_k = Q_tilde_X - Q_tilde_Y
U_n_k = apply(V_s_k, 2, sum)/sqrt(nrow(sample))
W_n_k = t(V_s_k)%*%V_s_k/nrow(sample)
}
T_n_k = as.numeric(t(U_n_k)%*%svd.inverse(W_n_k)%*%U_n_k)
statistics[j] <- T_n_k
}
k_optimal <- which.max(statistics - (1:max_j) * log(nrow(sample)))
T_tilde <- statistics[k_optimal]
return(T_tilde)
}

##### Quessy-Ethier (2012) test statistic
##### and its Bayesian bootstrap counterpart

QE.test = function(X,Y,n){

L = matrix(0,n,n)
for(p in 1:n){
 for(q in 1:n){
  L[p,q] = 2*max(X[p],Y[q]) - max(X[p],X[q]) - max(Y[p],Y[q])
 }
}
vec = seq(1,1,length = n)
QE = (vec %*% L %*% vec/n)[1,1]
return(QE)

}

QE.test.boot = function(X,Y,n,Xi){

L = matrix(0,n,n)
for(p in 1:n){
 for(q in 1:n){
  L[p,q] = 2*max(X[p],Y[q]) - max(X[p],X[q]) - max(Y[p],Y[q])
 }
}

vec.Xi = Xi/mean(Xi) - 1
QE.boot = (vec.Xi %*% L %*% vec.Xi/n)[1,1]
return(QE.boot)

}

##### Vexler et al. (2013) test statistic

VGH.test = function(X,Y,n,delta = 0.1){

Z = sort(X-Y)
alpha_n = n^(0.5+delta) 
beta_n = min( n^(1-delta), n/2 )
m = ceiling(alpha_n):floor(beta_n)
dm = length(m)
score = numeric( m[dm]-m[1]+1 )

for(i in 1:dm){
 V_mn = 1
  for(j in 1:n){ 
   d = 0
    for (k in 1:n){
      d = d + ifelse(Z[k] <= Z[ifelse(j + m[i] > n, n, j + m[i])], 1, 0) + 
        ifelse(-Z[k] <= Z[ifelse(j + m[i] > n, n, j + m[i])], 1, 0) - 
        ifelse(Z[k] <= Z[ifelse(j - m[i] < 1, 1, j - m[i])], 1, 0) - 
        ifelse(-Z[k] <= Z[ifelse(j - m[i] < 1, 1, j - m[i])], 1, 0)
    }
   Delta_jm <- d/(2*n)
   V_mn = V_mn*2*m[i]*(1-(m[i]+1)/(2*n))/(ifelse(Delta_jm == 0,1/n,Delta_jm)*n) 
  }
 score[i] <- V_mn
}

outcome = log(min(score))
return(outcome)

}

##### T test statistic (Konietschke and Pauly, 2014)

T.test = function(X,Y,n){

D = X-Y
T = sqrt(n)*ifelse(mean(D) == 0, 0, mean(D)/sd(D))
return(T)

}

##### New test statistics and their wild bootstrap counterparts

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

print(date())

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

valStat = matrix(0,9,Nk)

 P_W   = matrix(0,1,Nk)  
 P_BMP = matrix(0,1,Nk)
 P_MC  = matrix(0,1,Nk)
 P_QE  = matrix(0,1,Nk)
 P_VGH = matrix(0,1,Nk)
 P_KP  = matrix(0,1,Nk)
 P_KS  = matrix(0,1,Nk) 
 P_M   = matrix(0,1,Nk) 

for(i in 1:Nk){

##################################
### Example 1, Generator begins 
##################################

mu_X = 0
Sigma_X = 1

mu_Y = 0.8		### this parameter changes and it occurs in the figures as $mu$
Sigma_Y = 1

rho = -0.8/(Sigma_X * Sigma_Y)

Z1 = rnorm(n)
Z2 = rnorm(n)

X = Sigma_X * Z1 + mu_X
Y = Sigma_Y * (rho * Z1 + sqrt(1-rho^2) * Z2) + mu_Y

#################################
### Example 1, Generator ends 
#################################

		valStat[1,i] = W.test(X,Y,n)^2
	
	Z = c(X,Y)
	ldf1 = nparLD(Z~time, subject = subject, description = FALSE)
	sum.ldf1 = summary(ldf1)

		P_BMP[1,i] = ifelse( ifelse(sum.ldf1$Wald.test[1] == 0, 1, sum.ldf1$Wald.test[3]) < alpha, 1, 0)

	X.mat = matrix(0,n,2)
	X.mat[,1] = X		
	X.mat[,2] = Y

		P_MC[1,i] = MC.test(X.mat,hh,B0,B1)

	GPRY.val = GPRY.test(data.frame(X, Y))

### an asymptotic approximation for GPRY 

if(GPRY.val <= log(n)){
p.v = (2*pnorm(sqrt(GPRY.val))-1)*(2*pnorm(sqrt(log(n)))-1)
}else{
 if(GPRY.val >= 2*log(n)){
 p.v = (2*pnorm(sqrt(GPRY.val))-1)*(2*pnorm(sqrt(log(n)))-1) + 2*pnorm(-sqrt(log(n)))
 }else{
 lk = (2*pnorm(sqrt(GPRY.val))-1)*(2*pnorm(sqrt(log(n)))-1)
 pk = (2*pnorm(sqrt(GPRY.val))-1)*(2*pnorm(sqrt(log(n)))-1) + 2*pnorm(-sqrt(log(n)))
 a.par = (pk - lk)/log(n)
 p.v = a.par*GPRY.val + (pk - a.par*2*log(n))
}
}

		valStat[4,i] = ifelse((1-p.v) < alpha, 1, 0)
		valStat[5,i] = QE.test(X,Y,n)
		valStat[6,i] = VGH.test(X,Y,n,delta = 0.1)
		valStat[7,i] = T.test(X,Y,n)

	score = M.tests(X,Y,n,epsilon = 0.0001)

		valStat[8,i] = score[1]
		valStat[9,i] = score[4]

val_M_boot = matrix(0,4,Nk.boot)

 for(j in 1:Nk.boot){

	Eta = rexp(n)
		val_M_boot[1,j] = QE.test.boot(X,Y,n,Eta)

	Xi = rnorm(n)

		D.star = Xi*(X-Y)
		val_M_boot[2,j] = sqrt(n)*ifelse(mean(D.star) == 0, 0, mean(D.star)/sd(D.star))

	outcome = M.tests.boot(X,Y,n,epsilon = 0.0001,Xi,par.c)

		val_M_boot[3,j] = outcome[1]
		val_M_boot[4,j] = outcome[4]

 }	

T.p1 = mean( val_M_boot[2,] <= valStat[7,i] )
p.val.T = min(2*T.p1, 2-2*T.p1)

P_KP[1,i] = ifelse( p.val.T < alpha, 1,0)

p.v.KS = mean( valStat[8,i] < val_M_boot[3,] )
P_KS[1,i] = ifelse( p.v.KS < alpha, 1,0)

p.v.M = mean( valStat[9,i] < val_M_boot[4,] )
P_M[1,i] = ifelse( p.v.M < alpha, 1,0)

P_QE[1,i] = mean( val_M_boot[1,] > valStat[5,i] )

}

P_W[1,] = ifelse( valStat[1,] > qchisq(1-alpha,1),1,0 )

P_VGH[1,]  = ifelse( valStat[6,] > cv.VGH,1,0 )

##### Power functions

Pf.tests = matrix(0,1,9)  #

 Pf.tests[1,1] = mean(P_W[1,])*100   
 Pf.tests[1,2] = mean(P_BMP[1,])*100  
 Pf.tests[1,3] = mean(P_MC[1,]<alpha)*100 
 Pf.tests[1,4] = mean(valStat[4,])*100 
 Pf.tests[1,5] = mean(P_QE[1,]<alpha)*100 
 Pf.tests[1,6] = mean(P_VGH[1,])*100  
 Pf.tests[1,7] = mean(P_KP[1,])*100
 Pf.tests[1,8] = mean(P_KS[1,])*100
 Pf.tests[1,9] = mean(P_M[1,])*100

#####

print("Number of Monte Carlo runs")
print(Nk)

print("Number of wild bootstrap runs")
print(Nk.boot)

print("Number of observations")
print(n)

print("c")
print(par.c)

print("Power functions")

# print("W BMP MC GPRY QE VGH KP KS M")
# print(Pf.tests[1,])

print(cbind(c("W", "BMP", "MC", "GPRY", "QE", "VGH", "KP", "KS", "M"), Pf.tests[1,]))

print(date())
