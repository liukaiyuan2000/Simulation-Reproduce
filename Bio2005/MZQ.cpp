// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// notice that the code line 1 and line 2 are not explanation
#include <RcppArmadillo.h>
#include<math.h>
#include <cmath> 
#include <algorithm>  // max

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec CWnm(mat X){
  int N = X.n_rows;
  int n = N-1;
  int m = X.n_cols;
  mat R1 = cor(X);
  double Tn;
  int df;
  double pvalue;
  vec out(2);
  
  if(m<=n){
    Tn = -log(det(R1)*1.0)*(n-(2.0*m+5)/6);
  }else{
    Tn = 0;
  }
  Tn = max(Tn,1.0/m/n/m/n); 
  df=m*(m-1)/2;
  // pvalue <- 1-pchisq(Tn, df=df)
  //   out=list(Tn=Tn,pvalue=pvalue,df=df)
  out(0) = Tn;
  out(1) = df;
  
  return out;
}

// [[Rcpp::export]]
vec CSchottest(mat X){
  int N = X.n_rows;
  int n = N-1;
  int m = X.n_cols;
  mat R1 = cor(X);
  R1.diag() = zeros(R1.n_cols); 
  double Tn,esd,pvalue;
  vec out(2);
  
  Tn = accu(R1%R1)/2.0-m*(m-1)/2.0/n;
  esd = sqrt(m*(m-1)*(n-1.0)/(n+2.0)/n/n);
  Tn = Tn*1.0/esd;
  pvalue = 2*(1-normcdf(abs(Tn)));
  out(0) = Tn;
  out(1) = pvalue;
  
  return out;
}