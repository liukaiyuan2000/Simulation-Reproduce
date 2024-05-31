#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]



//' @title Calculate Tks.rand (rcpp version)
//'
//' @param V cbind(X.hat, Z) a n*(p+q) matrix
//' @param theta
//' @param e e.hat
//'
//' @return Tks test statistics(3.1)
//' @export
//'
//' @examples Tks_rand_rcpp(V, theta.normal, e.hat)
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
