#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::VectorXd Tks_rand_rcpp(Eigen::MatrixXd V, Eigen::MatrixXd theta, Eigen::MatrixXd e) {
  int n = e.rows();
  int rho = e.cols();
  int M = theta.cols();
  
  Eigen::MatrixXd Xt = V * theta;
  
  Eigen::MatrixXd Bn(n*M, rho);
  for (int i = 0; i < n*M; i++) {
    Eigen::MatrixXd diff = (Xt.array() <= Xt(i)).cast<double>();
    Eigen::MatrixXd prod = (e.transpose()*diff).array().abs();
    Bn.row(i) = prod.rowwise().sum();
  }
  
  Eigen::VectorXd Tn = Bn.colwise().maxCoeff() / (M * sqrt(n));
  
  return Tn;
}
