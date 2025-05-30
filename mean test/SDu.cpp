#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]

// Jiang Du
// May 22 2023
// Beijing University of Technology
// A test for the mean vector with fewer observations than the dimension. JMVA  
// Authors: SRIVASTAVA, M. S. and DU, M. (2008).


// [[Rcpp::export]]
double SDu(mat X,mat Y) {
  int n1=X.n_rows;
  int n2=Y.n_rows;
  int p=X.n_cols; 
  mat A1=cov(X);
  mat A2=cov(Y);
  mat Sn= ((n1-1.0)*A1+(n2-1.0)*A2)/(n1+n2-2.0);
  mat mX=mean(X,0);
  mat mY=mean(Y,0); 
  mat mu=mX-mY; 
  double n=n1+n2-2.0; 
  mat D2=zeros<mat>(p,p);
  D2.diag()=1.0/sqrt(Sn.diag());
  mat R=D2*Sn*D2;
  double cpn=1.0+trace(R*R)/((double) pow(p,3.0/2.0));
  double TN=0;
  for(int i=0;i<p;i++ ) TN=TN+pow(mu(i),2.0)/Sn(i,i);
  TN=1.0*(n1*n2)*TN/(n+2.0)- (1.0*n*p)/(n-2.0); 
  double Tn= (TN/sqrt(2.0*(trace(R*R)-1.0*pow(p,2)/n)*cpn)); 
  return(1.0-normcdf(Tn));
} 



