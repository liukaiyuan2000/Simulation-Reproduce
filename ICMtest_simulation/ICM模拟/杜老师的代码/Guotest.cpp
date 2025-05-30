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
  double cn = log(double(n))/double(n);
  vec lambda1 =arma::pow(lambda,2)+ cn;
  vec lambda2=zeros<vec>(p-1);
  for(int i=0;i<(p-1);i++){
    lambda2(i)= lambda1(i+1)/lambda1(i);
  }
  int q = lambda2.index_min();
  mat ALPHA(p,q+1);
  if(q==0) ALPHA=Alpha.col(0);
  if(q>=0){
    for(int i=0;i<=q;i++){
      ALPHA.col(i)=Alpha.col(i);
    }
  }
  return( ALPHA); 
}

// [[Rcpp::export]]
double Guotest(vec Y,mat X){
  mat Bn=CSE(Y,X);
  int n=Y.size(),p = X.n_cols;
  vec beta_hat=inv(X.t()*X)*X.t()*Y;
  vec res0=Y-X*beta_hat; 
  vec yhat0=X*beta_hat;
  mat index1=X*Bn; 
 double d=Bn.n_cols;
 double h=1.5*pow(double(n),-1.0/(4+d));
  double sigman=mean(res0%res0);
  
  mat W(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<=i;j++){
      W(i,j)=normcdf(norm((index1.row(i)-index1.row(j))/h));
      W(j,i)=W(i,j);
    }
  }
  W.diag() = zeros<vec>(n);
  double Tn=0,Tn1=0,Tn2=0;
  Tn1=sum(res0.t()*W*res0);
  vec res2=res0%res0;
  Tn2=2*sum(res2.t()*(pow(W,2))*res2);
  
  //cout<<Tn2<<endl;
  
  Tn= Tn1/sqrt(Tn2);
  return(Tn );
}






 


