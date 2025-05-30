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




// [[Rcpp::export]]
mat CSE(vec Y,mat X){
  int n = X.n_rows, p = X.n_cols;
  mat mu(p,1); 
  mu=(mean(X,0)).t();
  mat M;
  M.zeros(p,p);
  for(int i=0;i<n;i++){
    mat alpha(p,1); 
    for(int j=0;j<n;j++){
      alpha+=((X.row(j)-mu.t())*(Y(j)<=Y(i))).t();
    }
    alpha=alpha/n;
    M=M+alpha*alpha.t();
  }
  M=(M+M.t())/2; 
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, M);
  vec lambda=real(eigval);
  mat Alpha=eigvec; 
  return(Alpha); 
}






// [[Rcpp::export]]
vec ACM(vec Y,mat X,int B=500){
  mat Bn=CSE(Y,X);
  int n=Y.size();
  vec beta_hat=inv(X.t()*X)*X.t()*Y;
  vec res0=Y-X*beta_hat; 
  vec sigma_hat=sum(pow(res0,2))/n;
  vec yhat0=X*beta_hat;
  mat index1=X*Bn;
  double Tn=0;
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      Tn=Tn+res0(i)*res0(j)*exp(-norm(index1.row(i)-index1.row(j))/2-pow(yhat0(i)-yhat0(j),2)/2);
    }
  }
  Tn=Tn/n;
  vec TN=zeros<vec>(2);
  TN(0)=Tn;
  vec v,Ystar; 
  double pvalue=0;
  for(int b=0;b<B;b++){
    v=Mammencpp(n);
    Ystar=yhat0+(v%res0); 
    vec beta_hat_b=inv(X.t()*X)*X.t()*Ystar;
    vec resb=Ystar-X*beta_hat_b;
    vec yhatb=X*beta_hat_b;
    mat Bnb=CSE(Ystar,X);
    mat indexb=X*Bnb;
    double Tnb=0;
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        Tnb=Tnb+resb(i)*resb(j)*exp(-norm(indexb.row(i)-indexb.row(j))/2-pow(yhatb(i)-yhatb(j),2)/2);
      }
    } 
    Tnb=Tnb/n;
    if(Tnb>Tn) pvalue=pvalue+1.0;
  }
  TN(1)=pvalue/B; 
  return(TN);
  
}





