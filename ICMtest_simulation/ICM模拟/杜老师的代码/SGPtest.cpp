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





mat Wei(mat X){
  int n = X.n_rows;
  double temp=0;
  mat A(n,n); 
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      temp=0;
      for(int l=0;l<n;l++){
        if(all(X.row(i)<=X.row(l))){
          if(all(X.row(j)<=X.row(l))){
            temp=temp+1.0;
          }
        }
      }
      A(i,j)=temp;
    }
  }
  return( A/double(n)); 
}






// [[Rcpp::export]]
vec SGPtest(vec Y,mat X,int B=500){
  mat An=Wei(X);
  
  int n=Y.size();
  vec beta_hat=inv(X.t()*X)*X.t()*Y;
  vec res0=Y-X*beta_hat;
  vec yhat0=X*beta_hat;
  double Tn=as_scalar(res0.t()*An*res0);
  
 
  
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
    tmp=as_scalar(resb.t()*An*resb);
    //TN(b+1)=tmp;
    if(tmp>Tn) pvalue=pvalue+1.0;
  }
  TN(1)=pvalue/double(B);
  //return(pvalue/B);
  return(TN);
  
}





