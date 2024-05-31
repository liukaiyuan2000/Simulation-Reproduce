#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::VectorXd T3_rcpp(Eigen::MatrixXd X) {
  int N = X.rows();
  int n = N - 1;
  int m = X.cols();
  Eigen::MatrixXd centered = X.rowwise() - X.colwise().mean();
  Eigen::MatrixXd covMat = (centered.adjoint() * centered) / n;
  Eigen::MatrixXd covMatdiag0 = covMat;
  covMatdiag0.diagonal().setZero();
  double s2sum = covMat.diagonal().array().square().sum();
  double s4sum = covMat.diagonal().array().square().square().sum();
  double temp1 = covMatdiag0.array().square().sum();
  double temp2 = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      temp2 += covMat(i, i) * covMat(j, j) / n;
    }
  }
  temp2 = temp2 - s2sum;
  double lambda3 = (n / (n - 1.0)) * (temp1 - 1.0 / n * temp2) / s2sum + 1;
  double a20 = n * s2sum / m / (n + 2);
  double a40 = s4sum / m;
  double T3 = ((n - 1.0) * (lambda3 - 1.0)) / (2.0 * sqrt(1.0 - a40 / (m * pow(a20, 2))));
  double pvalue = 2 * (1 - R::pnorm(std::abs(T3), 0.0, 1.0, true, false));
  Eigen::VectorXd out(2);
  out << T3, pvalue;
  return out;
}
