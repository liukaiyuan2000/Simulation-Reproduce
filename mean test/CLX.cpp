#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]

// Jiang Du
// May 22 2024
// Beijing University of Technology
// Two-sample test of high dimensional means under dependence.  
// Authors: CAI, T. T., LIU, W. & XIA, Y. (2014).


// [[Rcpp::export]]
double CLX(mat X,mat Y) {
  int n1=X.n_rows;
  int n2=Y.n_rows;
  int p=X.n_cols; 
  mat A1=cov(X);
  mat A2=cov(Y);
  mat Sn= ((n1-1.0)*A1+(n2-1.0)*A2)/(n1+n2-2.0);
  mat mX=mean(X,0);
  mat mY=mean(Y,0); 
  mat mu=mX-mY; 
  vec TN = linspace(0, p-1, p);
  for(int i=0;i<p;i++ ) TN(i)=pow(mu(i),2.0)/Sn(i,i); 
  int n=n1+n2;
  double Tn=1.0*(n1*n2)*max(TN) /(1.0*n)-2.0*log(1.0*p)+log(log(1.0*p)); 
  double pv=1-exp(-exp(-Tn/2.0)/sqrt(3.1415926));
  return(pv);
} 



