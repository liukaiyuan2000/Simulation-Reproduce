#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::plugins(openmp)]]
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec Zhengtest(vec Y,mat X,double h){
  int n=Y.size(),p = X.n_cols;
  vec beta_hat=inv(X.t()*X)*X.t()*Y;
  vec res0=Y-X*beta_hat; 

  mat W(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<=i;j++){
      W(i,j)=normpdf(norm((X.row(i)-X.row(j))/h));
      W(j,i)=W(i,j);
    }
  }
  W.diag() = zeros<vec>(n);
  double Tn=0,Tn1=0,Tn2=0;
  Tn1=sum(res0.t()*W*res0);
  vec res2=res0%res0;
  Tn2=2*sum(res2.t()*(pow(W,2))*res2);
  //cout<<W<<endl;
//cout<<Tn1<<endl;
  //cout<<Tn2<<endl;
  Tn= Tn1/sqrt(Tn2);
  vec TN=zeros(2);
  TN(0)=Tn;
  TN(1)=2.0*(1-normcdf(double(abs(Tn))));
  return(TN);
}





