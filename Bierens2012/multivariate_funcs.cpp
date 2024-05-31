#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double Tnsc_rcpp(Eigen::MatrixXd x, Eigen::MatrixXd y, Eigen::MatrixXd ytilde, double c) {
  int m = y.cols();
  int k = x.cols();
  int n = x.rows();
  double I1 = 0;
  for(int j1 = 0; j1 < n - 1; j1 ++) {
    for(int j2 = j1 + 1; j2 < n; j2 ++) {
      double prody = 1;
      double prodytilde = 1;
      double prodyytilde = 1;
      double prodytildey = 1;
      double prodx = 1;
      for(int i = 0; i < m; i ++) {
        double diffy = c * (y(j1, i) - y(j2, i));
        double diffytilde = c * (ytilde(j1, i) - ytilde(j2, i));
        double diffyytilde = c * (y(j1, i) - ytilde(j2, i));
        double diffytildey = c * (ytilde(j1, i) - y(j2, i));
        prody = prody * (sin(diffy)/diffy);
        prodytilde = prodytilde * (sin(diffytilde)/diffytilde);
        prodyytilde = prodyytilde * (sin(diffyytilde)/diffyytilde);
        prodytildey = prodytildey * (sin(diffytildey)/diffytildey);
      }
      for(int i = 0; i < k; i ++) {
        double diffx = c * (x(j1, i) - x(j2, i));
        prodx = prodx * (sin(diffx)/diffx);
      }
      I1 += prodx * (prody + prodytilde - prodyytilde - prodytildey);
    }
  }
  double I2 = n;
  for(int j = 0; j < m; j ++) {
    double prod = 1;
    for(int i = 0; i < m; i ++) {
      double diff = c * (y(j, i) - ytilde(j, i));
      prod = prod * (sin(diff)/diff);
    }
    I2 -= prod;
  }
  double Tn = (2.0/n) * (I1 + I2);
  
  return(Tn);
}


// [[Rcpp::export]]
Eigen::VectorXd Tn_b_rcpp(
    Eigen::VectorXd x, Eigen::MatrixXd y_b, double c, int Model, int B
) {
  Environment myEnv = Environment::global_env();
  Function para_solve = myEnv["para.solve"];
  Function y_tilde_sample = myEnv["y.tilde.sample"];
  
  Eigen::VectorXd Tn(B);
  for(int i = 0; i < B; i ++) {
    P.middleCols<cols>(j)
    Eigen::MatrixXd y = y_b.middleCols<m>((i-1) * m);
    Eigen::MatrixXd ytilde = ytilde_b.middleCols<m>((i-1) * m);
    Tn(i) = Tnsc(x, y, ytilde, c);
  }
  
  return(Tn);
}







