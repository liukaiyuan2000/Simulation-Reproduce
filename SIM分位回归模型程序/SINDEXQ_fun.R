#main program for simulations
#contains: library, utility, data, estimation
#July 5, 2006; by Tracy Wu

##1.  invoke libraries

################################################################################
library(stats);
library(quantreg); #will call function lprq(), rq()
#library(lokern); # call function glkerns() to compute h_mean
#library(mda) #call for polyreg(), this is not the right polyreg; use own
library(KernSmooth);
#################################################################################

##2.  utility functions

#################################################################################
# 1
#  lprq0: local polynomial regression quantile. 
#          this estimates the quantile and its derivative at the point x_0
#  input: x - covariate sequence; y - response sequence; they form a sample.
#         h - bandwidth(scalar); tau - left-tail probability
#         x0 - point at which the quantile is estimated 
#  output: x0 - a scalar
#          fv - quantile est; dv - quantile derivative est
#################################################################################

lprq0<-function (x, y, h, tau = 0.5, x0)  #used in step 1 of the algorithm
{
    require(quantreg)
    fv <- x0
    dv <- x0
    
    z <- x - x0
    wx <- dnorm(z/h)
    r <- rq(y ~ z, weights = wx, tau = tau, ci = FALSE)
    fv <- r$coef[1]
    dv <- r$coef[2]
    list(x0 = x0, fv = fv, dv = dv)
}

##################################################################################
# 2
#  index.gamma(y,xx,p,gamma0,maxiter,crit)   #contains step 1 and 2
#  input: y - response vector; xx - covariate matrix; 
#         p - left-tail probability (quantile index), scalar
#         gamma0 - starting value of gamma 
#         maxiter - 
#         crit - 
#  output: gamma - index direction, unit norm, first nonneg component 
#          flag.conv  - whether the iterations converge
#  august 05
#### in step 2 of the proposed algorithm, we may consider random sampling a "subsample" of size, 
#      say 5n, of the augmented sample(with sample size n^2); 
#          do it several times and take average of the estimates 
#             (need to write a sampling step, but not in this utility function
####################################################################################
index.gamma<-function (y, xx, p, gamma0,maxiter,crit) 
{
flag.conv<-0; #flag whether maximum iteration is achieved

gamma.new<-gamma0; #starting value
gamma.new<-sign(gamma.new[1])*gamma.new/sqrt(sum(gamma.new^2));#单位化

n<-NROW(y); d<-NCOL(xx); 
a<-rep(0,n); b<-rep(0,n); #h<-rep(0,n);

iter<-1;
gamma.old<-2*gamma.new;
#gamma.old<-sign(gamma.old[1])*gamma.old/sqrt(sum(gamma.old^2));

while((iter < maxiter) & (sum((gamma.new-gamma.old)^2)>crit))
#while(iter < maxiter)
 {
 gamma.old<-gamma.new;
 iter<-iter+1;
 ####################################
 #  step 1: compute a_j,b_j; j=1:n  #
 ####################################
  a<-rep(0,n); b<-rep(0,n);#h<-rep(0,n);
  x<-rep(0,n);  
     for(jj in 1:d)
     {x<-x+xx[,jj]*gamma.old[jj]; #n-sequence, dim=null
       #x0<-x0+xx[j,jj]*gamma.old[jj]; #scalar
       } 
   hm<-dpill(x, y); #选择带宽
   hp<-hm*(p*(1-p)/(dnorm(qnorm(p)))^2)^.2;
 x0<-0;
 for(j in 1:n)
  {  
     x0<-x[j];
     fit<-lprq0(x, y, hp, p, x0) 
     a[j]<-fit$fv;
     b[j]<-fit$dv;
   }
 
 ############################# 
 # step 2: compute gamma.new #
 #############################
 # here, let v_j=1/n;
 ynew<-rep(0,n^2);
 xnew<-rep(0,n^2*d);
 xnew<-matrix(xnew,ncol=d);
 
 for (i in 1:n)
  { for (j in 1:n)
    { ynew[(i-1)*n+j]<-y[i]-a[j];
      for(jj in 1:d){ xnew[(i-1)*n+j,jj]<-b[j]*(xx[i,jj]-xx[j,jj]);}
    } 
  } 
 
  xg<-rep(0,n^2); #x*gamma
  for(jj in 1:d)
      {xg<-xg+xnew[,jj]*gamma.old[jj]; #n-sequence, dim=null
       }  
  xgh<-rep(0,n^2); #x*gamma/h
  for (i in 1:n)
  {   for (j in 1:n)
      { 
        xgh[(i-1)*n+j]<-xg[(i-1)*n+j]/hp;   
      } 
   } 
  wts<-dnorm(xgh);
  #fit<-rq(ynew ~0+ xnew, weights = wts, tau = p, method="fn") ; #pfn for very large problems
  fit<-rq(ynew ~0+ xnew, weights = wts, tau = p, ci = FALSE) ; #default: br method, for several thousand obs
          # 0, to exclude intercept
  gamma.new<-fit$coef;
  gamma.new<-sign(gamma.new[1])*gamma.new/sqrt(sum(gamma.new^2));   #normalize

} #end iterates over iter;

iter<-iter;

flag.conv<-(iter < maxiter) ;# = 1 if converge; =0 if not converge
#flag.conv<- 1- ((iter=maxiter)&(sum((gamma.new-gamma.old)^2)<crit))

gamma<-gamma.new;
list(gamma=gamma,flag.conv=flag.conv)
}


##################################################################################
# 3 
#  index.gamma1(y,xx,p,gamma0,maxiter,crit,subsize,subrep)   #contains step 1 and 2
#     it's revised version of index.gamma(), in that random sampling is implemented in its step 2
#
#  input: y - response vector; xx - covariate matrix; 
#         p - left-tail probability (quantile index), scalar
#         gamma0 - starting value of gamma 
#         maxiter - 
#         crit - 
#  output: gamma - index direction, unit norm, first nonneg component 
#          flag.conv  - whether the iterations converge
#  august 05
#
### in step 2 of the proposed algorithm, we consider random sampling a "subsample" of size, 
#     say 5n, of the augmented sample(with sample size n^2); 
#          do it several times and take average of the estimates 
#            
####################################################################################

index.gamma1<-function (y, xx, p, gamma0,maxiter,crit,subsize,subrep) 
{
flag.conv<-0; #flag whether maximum iteration is achieved

gamma.new<-gamma0; #starting value
gamma.new<-sign(gamma.new[1])*gamma.new/sqrt(sum(gamma.new^2));

n<-NROW(y); d<-NCOL(xx); 
a<-rep(0,n); b<-rep(0,n); #h<-rep(0,n);
gamma.inter<-rep(0,d*subrep); 

iter<-1;
gamma.old<-2*gamma.new;
#gamma.old<-sign(gamma.old[1])*gamma.old/sqrt(sum(gamma.old^2));

while((iter < maxiter) & (sum((gamma.new-gamma.old)^2)>crit))
#while(iter < maxiter)
 {
 gamma.old<-gamma.new;
 iter<-iter+1;
 ####################################
 #  step 1: compute a_j,b_j; j=1:n  #
 ####################################
  a<-rep(0,n); b<-rep(0,n);#h<-rep(0,n);
  x<-rep(0,n);  
     for(jj in 1:d)
     {x<-x+xx[,jj]*gamma.old[jj]; #n-sequence, dim=null
       #x0<-x0+xx[j,jj]*gamma.old[jj]; #scalar
       } 
   hm<-dpill(x, y); 
   hp<-hm*(p*(1-p)/(dnorm(qnorm(p)))^2)^.2;
 x0<-0;
 for(j in 1:n)
  {  
     x0<-x[j];
     fit<-lprq0(x, y, hp, p, x0) 
     a[j]<-fit$fv;
     b[j]<-fit$dv;
   }
 
 ############################# 
 # step 2: compute gamma.new #
 #############################
 # here, let v_j=1/n;
 # augmented "observations"
 ynew<-rep(0,n^2); 
 xnew<-rep(0,n^2*d);
 xnew<-matrix(xnew,ncol=d);
 for (i in 1:n)
    { for (j in 1:n)
      { ynew[(i-1)*n+j]<-y[i]-a[j];
      for(jj in 1:d){ xnew[(i-1)*n+j,jj]<-b[j]*(xx[i,jj]-xx[j,jj]); }
     } 
   }  #generated ynew and xnew
  xg<-rep(0,n^2); #x*gamma
  for(jj in 1:d)
      {xg<-xg+xnew[,jj]*gamma.old[jj]; #n-sequence, dim=null
       }  
  xgh<-rep(0,n^2); #x*gamma/h
  for (i in 1:n)
  {   for (j in 1:n)
      { 
        xgh[(i-1)*n+j]<-xg[(i-1)*n+j]/hp;   
      } 
   } 
 wts<-dnorm(xgh); #generated wts, the weights
   ######################
   #  sampling step
   ######################
 gamma.inter<-rep(0,d*subrep); 
 gamma.inter<-matrix(gamma.inter, ncol=subrep);  #to store intermediate gammas
 for (subiter in 1:subrep)
  {   id.sample<-  sample(seq(1,n^2), subsize) #default: w/o replacement  
      subx<-xnew[id.sample,]
      suby<-ynew[id.sample]
      subw<-wts[id.sample]

      subfit<-rq(suby ~0+ subx, weights = subw, tau = p, ci = FALSE) ; #default: br method, for several thousand obs
          # 0, to exclude intercept
    subgamma<-subfit$coef;
    subgamma<-sign(subgamma)*subgamma/sqrt(sum(subgamma^2));   #normalize
    gamma.inter[,subiter]<-subgamma;
  }  #end iterates over subiter
 ###################################
 #take average of gamma estimates
 ###################################
 gamma.inter<-gamma.inter #dim: d by subrep
 gamma.new<-apply(gamma.inter,1,mean)
 
 gamma.new<- sign(gamma.new)*gamma.new/sqrt(sum(gamma.new^2))  

} #end iterates over iter;

gamma.new<-gamma.new
iter<-iter;

flag.conv<-(iter < maxiter) ;# = 1 if converge; =0 if not converge
#flag.conv<- 1- ((iter=maxiter)&(sum((gamma.new-gamma.old)^2)<crit))

gamma<-gamma.new;
list(gamma=gamma,flag.conv=flag.conv)
}

####################
#end of index.gamma1
#####################




  
    

