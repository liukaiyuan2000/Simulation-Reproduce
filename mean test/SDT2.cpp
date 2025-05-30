#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
vec SDT2(mat X,mat Y) {
  int n1=X.n_rows;
  int n2=Y.n_rows;
  int p=X.n_cols; 
  mat A1=cov(X);
  mat A2=cov(Y);
  mat Sn= ((n1-1.0)*A1+(n2-1.0)*A2)/(n1+n2-2.0);
  mat mX=mean(X,0);
  mat mY=mean(Y,0); 
  mat mu=mX-mY; 
 double n=1.0*n1+1.0*n2;
 double a1=trace(Sn)/double(p);
 double a2=pow(n,2.0)*((trace(Sn*Sn))-pow(trace(Sn),2.0)/(1.0*n))/(1.0*n-1.0)/(1.0*n+2.0)/double(p);
 double rhat=double(p)*pow(a1,2.0)/a2;
 double Tn=1.0*(n1*n2)*sum(sum(pow(mu,2.0)))/(1.0*n)/trace(Sn); 
 vec TN=zeros<vec>(3);
 TN(0)=Tn;
 TN(1)= floor(rhat);
 TN(2)= floor(n*rhat);
  return(TN);
} 



