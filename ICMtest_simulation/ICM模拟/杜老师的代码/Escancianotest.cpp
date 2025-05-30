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



double angle_xy(rowvec x, rowvec y) {
  double xnorm = norm(x);
  double ynorm = norm(y);
  if(xnorm == 0 || ynorm == 0) {
    return 0;
  } else {
    double cxy = norm_dot(x, y);
    if(cxy > 1) cxy = 1;
    if(cxy < -1) cxy = -1;
    return acos(cxy);
  }
}

double pconstant(int p) {
  const double pi = 3.14159265358979323846;
  double l1 = pow(pi, p/2.0 - 1);
  double l2 = ::tgamma(p/2.0+1);
  return l1 / l2;
}

//[[Rcpp::export]]
mat angle_matrix(mat x) {
  const double pi = 3.14159265358979323846;
  int n = x.n_rows, p = x.n_cols;
  double pcon = pconstant(p + 1);
  cube am(n, n, n);
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      for (int j = i; j < n; ++j) {
        am(i, j, k) = pcon * (pi - angle_xy(x.row(i) - x.row(k), x.row(j) - x.row(k)));
        if(i != j) am(j, i, k) = am(i, j, k);
      }
    }
  }
  
  mat angm = sum(am, 2);
  return angm;
}





//[[Rcpp::export]]
mat angle_matrix_parallel(mat x, int ncores = 3) {
  const double pi = 3.14159265358979323846;
  int n = x.n_rows, p = x.n_cols;
  double pcon = pconstant(p + 1);
  cube am(n, n, n);
#pragma omp parallel for num_threads(ncores)
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      for (int j = i; j < n; ++j) {
        am(i, j, k) = pcon * (pi - angle_xy(x.row(i) - x.row(k), x.row(j) - x.row(k)));
        if(i != j) am(j, i, k) = am(i, j, k);
      }
    }
  }
  
  mat angm = sum(am, 2);
  return angm;
}




vec indexUvec(vec x,uvec inds){
  int n=x.size();
  vec newVec = zeros<vec>(n);
  for (int i = 0; i < n; i++){
    newVec(i) = x(inds(i));
  }
  return newVec;
}



 

// [[Rcpp::export]]
vec Escanciano(vec Y,mat X,int B=500){
  int n=Y.size();
  mat A=angle_matrix_parallel(X);
  
  vec beta_hat=inv(X.t()*X)*X.t()*Y;
  vec res0=Y-X*beta_hat; 
  vec yhat0=X*beta_hat;
  double Tn=as_scalar(res0.t()*A*res0/n/n/n);
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
    tmp=as_scalar(resb.t()*A*resb/n/n/n);
    //TN(b+1)=tmp;
    if(tmp>Tn) pvalue=pvalue+1.0;
  }
  TN(1)=pvalue/B;
  //return(pvalue/B);
  return(TN);
  
}



