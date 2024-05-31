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

Eigen::MatrixXd scor_rcpp(Eigen::MatrixXd X) {
  Function base_cor("cor");
  NumericMatrix res(base_cor(X, Named("method") = "spearman"));
  Eigen::MatrixXd result(as<Eigen::MatrixXd>(res));
  
  return result;
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
  double trs2 = (covMat* covMat).trace();
  double trs = covMat.trace();
  double lambda3 = (n / (n - 1.0)) * (trs2 - pow(trs, 2) / n) / s2sum;
  double a201 = n * s2sum / m / (n + 2.0);
  double a202 = s2sum / m;
  double a40 = s4sum / m;
  double v1 = 1.0 - a40 / (m * pow(a201, 2));
  double v2 = 1.0 - a40 / (m * pow(a202, 2));
  double v;
  if(v1 < 0){
    v = v2;
  }else{
    v = v1;
  };
  double T3 = (n * (lambda3 - 1.0)) / (2.0 * sqrt(v));
  double pvalue = 2 * (1 - R::pnorm(std::abs(T3), 0.0, 1.0, true, false));
  Eigen::VectorXd out(2);
  out << T3, pvalue;
  return out;
}


// [[Rcpp::export]]
Eigen::VectorXd Trho_rcpp(Eigen::MatrixXd X) {
  int N = X.rows();
  int n = N - 1;
  int m = X.cols();
  Eigen::MatrixXd R1 = scor_rcpp(X);
  R1.diagonal().setZero();
  
  Eigen::VectorXd out(2);
  double Trho = R1.array().square().sum() / 2.0 - m * (m - 1) / 2.0 / n;
  double esd = m / N;
  Trho = Trho / esd;
  double pvalue = 2 * (1 - R::pnorm(std::abs(Trho), 0.0, 1.0, true, false));
  
  out << Trho, pvalue;
  return out;
}


// [[Rcpp::export]]
Eigen::VectorXd Snm_rcpp(Eigen::MatrixXd X) {
  int N = X.rows();
  int n = N - 1;
  int m = X.cols();
  Eigen::MatrixXd R1 = scor_rcpp(X);
  R1.diagonal().setZero();
  
  Eigen::VectorXd out(2);
  double Sn = R1.array().square().sum() / 2.0 - m * (m - 1) / 2.0 / n;
  double esd = sqrt(m * (m - 1.0) * (25 * pow(N, 3) - 57 * pow(N, 2) - 40 * N + 108) / (25 * pow(n, 3) * N * (N + 1)));
  Sn = Sn / esd;
  double pvalue = 2 * (1 - R::pnorm(std::abs(Sn), 0.0, 1.0, true, false));
  
  out << Sn, pvalue;
  return out;
}


