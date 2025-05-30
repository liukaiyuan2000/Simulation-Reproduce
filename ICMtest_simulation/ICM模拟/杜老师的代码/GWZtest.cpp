#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::plugins(openmp)]]
//[[Rcpp::depends(RcppArmadillo)]]





// [[Rcpp::export]]
mat CSEcpp1(vec Y,mat X){
  int n = X.n_rows, p = X.n_cols;
  mat A=cov(X); 
  mat mu(p,1);
  mu=(mean(X,0)).t(); 
  mat M=zeros(p,p);
  mat alpha(p,1); 
  double a1=0;
  for(int i=0;i<n;i++){
    alpha=zeros(p,1); 
    for(int j=0;j<n;j++){
      a1=0;
      if(Y(j)<=Y(i)) a1=1.0;
      alpha+=a1*((X.row(j)-mu.t())).t();
    }
    alpha=alpha/double(n);
    M=M+alpha*alpha.t();
  }
  
  M=M/double(n);
  cx_vec eigval;
  cx_mat eigvec; 
  vec  lambda;
  eig_gen(eigval, eigvec, inv(A)*M);
  lambda=real(eigval);
  mat Alpha=zeros(p,p);
  double cn =2.0*pow(double(n),3.0/4.0)/double(p);
  vec lambda1=lambda%lambda;
  
  lambda1=double(n)*cumsum(lambda1)/sum(lambda1)-cn*linspace(1, p, p)%linspace(2, p+1, p)/2.0;
  int q = lambda1.index_max();
  mat ALPHA(p,q+1);
  if(q==0) ALPHA=real( eigvec.col(0));
  if(q>=0){
    for(int i=0;i<=q;i++){
      ALPHA.col(i)=real( eigvec.col(i));
    }
  }
  return( ALPHA); 
}





// [[Rcpp::export]]
double GWZtest(vec Y,mat X){
  mat Bn=CSEcpp1(Y,X);
  
  int n=Y.size();
  vec beta_hat=inv(X.t()*X)*X.t()*Y;
  vec res0=Y-X*beta_hat; 
  mat index1=X*Bn; 
  double q=Bn.n_cols;
  double h=1.5*pow(double(n),-1.0/(4+q));
  mat W(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<=i;j++){
      W(i,j)=normpdf(norm(index1.row(i)-index1.row(j))/h);
      W(j,i)=W(i,j);
    }
  }
  W.diag() = zeros<vec>(n);
  double Tn=0,Tn1=0,Tn2=0;
  Tn1=sum(res0.t()*W*res0);
  vec res2=res0%res0;
  Tn2=2.0*sum(res2.t()*(pow(W,2))*res2);
  
  Tn= Tn1/sqrt(Tn2);
  return(Tn );
}









