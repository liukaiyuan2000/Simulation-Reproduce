#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::plugins(openmp)]]
//[[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
vec Mammencpp(int n){
  srand((unsigned)time(NULL)); 
  double  p=(sqrt(5.0)+1.0)/(2.0*sqrt(5.0));
  vec  v = zeros<vec>(n)-(sqrt(5.0)-1.0)/2.0;
  
  vec  b = randu<vec>(n);
  
  for(int i=0;i<n;i++){
    if(b(i)>p) v(i)=(sqrt(5)+1)/2;
  }
  return(v); 
}



 

//[[Rcpp::export]]
mat ICMA(mat x) {
  int n = x.n_rows;
  mat A(n,n);
  for (int i = 0; i < n; ++i) {
      for (int j = i; j < n; ++j) {
        A(i, j) = norm(x.row(i) - x.row(j));
        A(j,i)=A(i,j);
      }
    }
  return A;
}


// [[Rcpp::export]]
vec ICM(vec Y,mat X,int B=500){
  int n=Y.size();
  mat A=ICMA(X);
  
  vec beta_hat=inv(X.t()*X)*X.t()*Y;
  vec res0=Y-X*beta_hat; 
  vec yhat0=X*beta_hat;
  mat AA=exp(-0.5*pow(A,2));
  double Tn=as_scalar(res0.t()*AA*res0/n);
  vec v,Ystar; 
  double pvalue=0;
  double tmp=0;
  vec TN=zeros<vec>(2);
  TN(0)=Tn;
  for(int b=0;b<B;b++){
    v=Mammencpp(n);
    Ystar=yhat0+(v%res0); 
    vec beta_hat_b=inv(X.t()*X)*X.t()*Ystar;
    vec resb=Ystar-X*beta_hat_b;
    //cout<<beta_hat_b<<endl;
    tmp=as_scalar(resb.t()*AA*resb/n);
    //TN(b+1)=tmp;
    if(tmp>Tn) pvalue=pvalue+1.0;
  }
  TN(1)=pvalue/B;
  //cout<<beta_hat<<endl;
  //cout<<A<<endl;
  //return(pvalue/B);
  return(TN);
  
}





 