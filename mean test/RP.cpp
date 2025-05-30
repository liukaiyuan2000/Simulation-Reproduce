#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]

// Jiang Du
// May 22 2023
// Beijing University of Technology
//  A more powerful two-sample test in high dimensions using random projection 
// Authors:LOPES, M. E., JACOB, L. and WAINWRIGHT, M. J.



// [[Rcpp::export]]
double RP(mat X,mat Y) {
  int n1=X.n_rows;
  int n2=Y.n_rows;
  int p=X.n_cols; 
  int n=n1+n2-2;
  int k=floor(n/2);
  mat Pk(p, k, fill::randn);
  mat A1=cov(X);
  mat A2=cov(Y);
  mat Sn= ((n1-1.0)*A1+(n2-1.0)*A2)/(n1+n2-2.0);
  mat P=Pk*inv(Pk.t()*Sn*Pk)*Pk.t();
  mat mX=mean(X,0);
  mat mY=mean(Y,0); 
  mat mu=mX-mY; 
  
   
  double Tn=as_scalar(1.0*(n1*n2)*mu*P*mu.t()/(n+2.0)); 
 
  return(Tn);
} 



