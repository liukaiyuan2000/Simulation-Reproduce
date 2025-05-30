#include <RcppArmadillo.h>
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]
 
// Jiang Du
// May 22 2023
// Beijing University of Technology
// A two-sample test for high-dimensional data with applications to gene-set testing. Ann. Statist. 38 808?C835.
// Authors: CHEN, S. X. and QIN, Y.-L.


double XCQ(mat X) {
 int n1=X.n_rows; 
 int p=X.n_cols; 
 mat Tn = zeros<mat>(p,p);
 mat Xjk=sum(X); 
 for(int j=0;j<n1;j++){
   for(int k=0;k<n1;k++){ 
     if(k!=j){
       Tn = Tn+ X.row(j).t()*(X.row(j)-(Xjk-X.row(j)-X.row(k))/(1.0*n1-2.0))*X.row(k).t()*(X.row(k)-(Xjk-X.row(j)-X.row(k))/(1.0*n1-2.0));  
     }
   }
 } 
 double num1=1.0*n1*(n1-1);
 return( trace(Tn)/num1);
}


double XYCQ(mat X,mat Y) {
 int n1=X.n_rows; 
 int n2=Y.n_rows; 
 int p=X.n_cols; 
 mat Tn = zeros<mat>(p,p);
 mat Xjk=sum(X); 
 mat Yjk=sum(Y); 
 for(int j=0;j<n1;j++){
   for(int k=0;k<n2;k++){  
     Tn = Tn+ X.row(j).t()*(X.row(j)-(Xjk-X.row(j))/(1.0*n1-1.0))*Y.row(k).t()*(Y.row(k)-(Yjk-Y.row(k))/(1.0*n1-1.0));  
     
   }
 } 
 double num1=1.0*n1*n2;
 return( trace(Tn)/num1);
}




// [[Rcpp::export]]
double CQ(mat X,mat Y) {
 int n1=X.n_rows;
 int n2=Y.n_rows;
  
 double  A1 =0;
 double  A2 =0;
 double  A3 =0;
 
 for(int i=1;i<n1;i++){
   for(int j=0;j<i;j++){
     A1 = A1+sum(X.row(i)%X.row(j));  
   }
 }
 for(int i=1;i<n2;i++){
   for(int j=0;j<i;j++){
     A2  = A2+sum(Y.row(i)%Y.row(j));  
   }
 }
 
 for(int i=1;i<n1;i++){
   for(int j=0;j<n2;j++){
     A3  = A3+sum(X.row(i)%Y.row(j));  
   }
 }
 
 double num1=1.0*n1*(n1-1);
 double num2=1.0*n2*(n2-1);
 double num3=1.0*n1*n2;
 double esd=sqrt(2*XCQ(X)/num1+2*XCQ(Y)/num2+4*XYCQ(X,Y)/num3);  
 double Tn=2.0*A1/(1.0*(pow(n1,2.0)-n1))+2.0*A2/(1.0*(pow(n2,2.0)-n2))-2*A3/(1.0*(n1*n2));
 return(1.0-normcdf(Tn/esd)); 
} 
 
 
 
 