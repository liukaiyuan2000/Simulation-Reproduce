#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::plugins(openmp)]]
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double SZtest(vec Y,mat X){
  int n=Y.size();
  vec beta_hat=inv(X.t()*X)*X.t()*Y;
  vec res0=Y-X*beta_hat; 
  double sigman=0;
  for(int i=0;i<n;i++){
    sigman=sigman+pow(res0(i),2);
  }
  sigman=sigman/double(n);
  vec index=X*beta_hat;
  vec a0=zeros<vec>(3);
  a0(0)=0.97;a0(1)=0.98;a0(2)=0.99;
  vec x10=quantile(sort(index),a0);
  double x0=x10(2);
  mat W(n,n),tW(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(index(i)<=index(j)){
        W(i,j)=1;
        tW(i,j)=0;
        }else{
          W(i,j)=0;
          tW(i,j)=1;
        }
      }
    }
  double sigma2=as_scalar(mean(pow(res0,2)));
  vec Rn1=zeros<vec>(n),An=zeros<vec>(n);
  vec  alpha=index/sigma2;
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      An(i)=An(i)+sigma2*pow(alpha(j),2)*W(i,j)/(double(n));
      Rn1(i)= Rn1(i)+W(j,i)*res0(j)/sqrt(double(n));
      }
    }
  vec Rn2=zeros<vec>(n),TRn=zeros<vec>(n); 
  mat temp1(n,1),temp2(n,1);
  temp2.col(0)=res0%index;
  vec sgn=zeros<vec>(n);
  for(int i=0;i<n;i++){
    temp1=W.col(i)%index/An;
    Rn2(i)=as_scalar(sum(temp1.t()*W*temp2)/sqrt(double(n))/double(n));
    TRn(i)=Rn1(i)-Rn2(i);
    if(index(i)<=x0) sgn(i)=1;
    }
  double Tn=0;
  Tn=mean(pow(TRn,2)%sgn )/mean(sgn)/mean(sgn)/sigma2;
  return(Tn );
}





