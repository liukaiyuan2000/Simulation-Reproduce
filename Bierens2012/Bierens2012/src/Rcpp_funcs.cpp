#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double Tnsc_rcpp(Eigen::VectorXd x, Eigen::VectorXd y, Eigen::VectorXd ytilde, double c) {
  int n = x.rows();
  double I1 = 0;
  for(int j1 = 0; j1 < n - 1; j1 ++) {
    for(int j2 = j1 + 1; j2 < n; j2 ++) {
      double diffy = c * (y(j1) - y(j2));
      double diffytilde = c * (ytilde(j1) - ytilde(j2));
      double diffyytilde = c * (y(j1) - ytilde(j2));
      double diffytildey = c * (ytilde(j1) - y(j2));
      double prody, prodytilde, prodyytilde, prodytildey;
      if (diffy == 0) {
        prody = 1;
      } else{
        prody = sin(diffy)/diffy;
      }
      if (diffytilde == 0) {
        prodytilde = 1;
      } else{
        prodytilde = sin(diffytilde)/diffytilde;
      }
      if (diffyytilde == 0) {
        prodyytilde = 1;
      } else{
        prodyytilde = sin(diffyytilde)/diffyytilde;
      }
      if (diffytildey == 0) {
        prodytildey = 1;
      } else{
        prodytildey = sin(diffytildey)/diffytildey;
      }
      double diffx = c * (x(j1) - x(j2));
      double prodx = sin(diffx)/diffx;

      I1 += prodx * (prody + prodytilde - prodyytilde - prodytildey);
    }
  }

  Eigen::VectorXd diff = c * (y - ytilde);
  Eigen::VectorXd temp1 = Eigen::VectorXd::Ones(n);
  for (int i = 0; i < n; i++) {
    if (diff(i) != 0) {
      temp1(i) = sin(diff(i)) / diff(i);
    }
  }
  // std::cout << temp1;
  double temp2 = temp1.array().sum();
  double I2 = n - temp2;
  double Tn = (2.0/n) * (I1 + I2);

  return(Tn);
}
