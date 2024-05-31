#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <random>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace Eigen;
using namespace std;

// [[Rcpp::export]]
NumericMatrix mvrnorm_rcpp(int n, NumericVector mu, NumericMatrix sigma) {
  int p = mu.size();
  NumericMatrix X(n, p);
  NumericMatrix L = clone(sigma);
  
  for (int i = 0; i < p; i++) {
    for (int j = 0; j <= i; j++) {
      for (int k = 0; k < j; k++) {
        L(i, j) -= L(i, k) * L(j, k);
      }
      if (i == j) {
        L(i, j) = sqrt(L(i, j));
      } else {
        L(i, j) /= L(j, j);
      }
    }
  }
  
  for (int i = 0; i < n; i++) {
    NumericVector Z = rnorm(p);
    for (int j = 0; j < p; j++) {
      X(i, j) = mu[j];
      for (int k = 0; k <= j; k++) {
        X(i, j) += L(j, k) * Z[k];
      }
    }
  }
  
  return X;
}

Eigen::MatrixXd cor_rcpp(Eigen::MatrixXd X) {
  int n = X.rows();
  Eigen::MatrixXd centered = X.rowwise() - X.colwise().mean();
  Eigen::MatrixXd covMat = (centered.adjoint() * centered) / (n - 1);
  Eigen::VectorXd sd = covMat.diagonal().cwiseSqrt();
  Eigen::MatrixXd corMat = covMat.array().colwise() / sd.array();
  corMat = corMat.array().rowwise() / sd.transpose().array();
  return corMat;
}


// [[Rcpp::export]]
Eigen::VectorXd Wnm_rcpp(Eigen::MatrixXd X) {
  int N = X.rows();
  int n = N - 1;
  int m = X.cols();
  Eigen::MatrixXd R1 = cor_rcpp(X);
  double Tn,pvalue;
  int df;
  Eigen::VectorXd out(2);
  
  if(m <= n) {
    Tn = -log(R1.determinant() * 1.0) * (n - (2.0 * m + 5.0) / 6.0);
  } else {
    Tn = 0;
  }
  Tn = std::max(Tn, 1.0 / (m * n * m * n));
  df = m * (m - 1) / 2;
  pvalue = 1 - R::pchisq(Tn, df, true, false);
  out << Tn, pvalue;
  return out;
}

// [[Rcpp::export]]
Eigen::VectorXd Tnm_rcpp(Eigen::MatrixXd X) {
  int N = X.rows();
  int n = N - 1;
  int m = X.cols();
  Eigen::MatrixXd R1 = cor_rcpp(X);
  R1.diagonal().setZero();
  
  Eigen::VectorXd out(2);
  double Tn = R1.array().square().sum() / 2.0 - m * (m - 1) / 2.0 / n;
  double esd = sqrt(m * (m - 1.0) * (n - 1) / (n + 2) / n / n);
  Tn = Tn / esd;
  double pvalue = 2 * (1 - R::pnorm(std::abs(Tn), 0.0, 1.0, true, false));
  
  out << Tn, pvalue;
  return out;
}





