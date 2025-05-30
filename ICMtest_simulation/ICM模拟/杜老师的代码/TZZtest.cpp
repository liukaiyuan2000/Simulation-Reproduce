#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::plugins(openmp)]]
//[[Rcpp::depends(RcppArmadillo)]]








// [[Rcpp::export]]
mat CSEcpp(vec Y,mat X){
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
vec BTn(int B=500,int M=500,double T1=1){
  
  vec TN=zeros(B);
  double dt=sqrt(T1/(double(M)-1.0));//To ensure we start at zero and end at 1
  for(int i=0;i<B;i++){
    vec W=zeros<vec>(M);
    vec N=randn(M);
    for(int j=1;j<M;j++){
      W(j)=W(j-1)+dt*N(j);
    }
    TN(i)=mean(W%W);
    
  }
  
  
  return(TN);
}





// [[Rcpp::export]]
vec TZZtest(vec Y,mat X){
  mat Beta1=CSEcpp(Y,X);
  vec beta1=Beta1.col(0);
  int B=500,M=500;
  vec WBN=BTn(B,M);
  int n=Y.size();
  vec beta_hat=inv(X.t()*X)*X.t()*Y;
  vec res0=Y-X*beta_hat; 
  
  vec index=X*beta1;
  //vec res0=Y-0.1821531-1.1419628*index+0.3577954*pow(index,2)+0.1638841*pow(index,3);
  vec a0=zeros<vec>(3);
  a0(0)=0.97;a0(1)=0.98;a0(2)=0.99;
  vec x10=quantile(sort(index),a0);
  double x0=x10(2);
  mat W(n,n),tW(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(index(i)<=index(j)){
        W(i,j)=1.0;
        tW(i,j)=0;
      }else{
        W(i,j)=0;
        tW(i,j)=1.0;
      }
    }
  }
  double sigma2=as_scalar(mean(res0%res0));
  vec Rn1=zeros<vec>(n),An=zeros<vec>(n);
  vec  alpha=index/sigma2;
  for(int i=0;i<n;i++){
    An(i)=sigma2*sum((alpha%alpha).t()%W.row(i))/(double(n));
    Rn1(i)=sum(W.col(i)%res0)/sqrt(double(n));
  }
  vec Rn2=zeros<vec>(n),TRn=zeros<vec>(n); 
  mat temp1(n,1),temp2(n,1);
  temp2.col(0)=res0%alpha;
  vec sgn=zeros<vec>(n);
  for(int i=0;i<n;i++){
    temp1=W.col(i)%index/An;
    Rn2(i)=as_scalar(sum(temp1.t()*W*temp2)/sqrt(double(n))/double(n));
    TRn(i)=Rn1(i)-Rn2(i);
    if(index(i)<=x0) sgn(i)=1.0;
  }
  double Tn=0;
  Tn=mean(pow(TRn,2.0)%sgn )/mean(sgn)/mean(sgn)/sigma2;
  vec TNN=zeros<vec>(2);
  TNN(0)=Tn;
  for(int i=0;i<B;i++){
    if(WBN(i)>=Tn){
      TNN(1)=TNN(1)+1.0;
    }
  }
  TNN(1)=TNN(1)/double(B);
  return(TNN );
}





