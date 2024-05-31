#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat aux_func(arma::vec x, arma::vec y, double c) {
  int n = x.n_elem;
  arma::mat result = arma::mat(n, n, arma::fill::ones) * M_PI;

  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      result(i, j) = (x(i) - y(j)) * c;
    }
  }

  return result;
}

//' @title A Rcpp version for Tnsc
//' @export
// [[Rcpp::export]]
double Tnsc_rcpparmadillo(arma::vec x, arma::vec y, arma::vec ytilde, double c) {
 int n = x.n_elem;
 arma::mat diffx = aux_func(x, x, c);
 arma::mat diffy = aux_func(y, y, c);
 arma::mat diffytilde = aux_func(ytilde, ytilde, c);
 arma::mat diffyytilde = aux_func(y, ytilde, c);
 arma::mat diffytildey = aux_func(ytilde, y, c);

 arma::mat I11 = arma::sin(diffy) / diffy;
 arma::mat I12 = arma::sin(diffytilde) / diffytilde;
 arma::mat I13 = arma::sin(diffyytilde) / diffyytilde;
 arma::mat I14 = arma::sin(diffytildey) / diffytildey;
 arma::mat I15 = arma::sin(diffx) / diffx;
 I11.replace(arma::datum::nan, 1);
 I12.replace(arma::datum::nan, 1);
 I13.replace(arma::datum::nan, 1);
 I14.replace(arma::datum::nan, 1);
 arma::mat I1 = (I11 + I12 - I13 - I14) % I15;
 double sumI1 = arma::accu(I1);
 arma::vec diff = c * (y - ytilde);
 arma::vec temp1 = arma::sin(diff) / diff;
 temp1.replace(arma::datum::nan, 1);
 double temp2 = arma::accu(temp1);
 double I2 = n - temp2;
 double Tn = (2.0/n) * (sumI1 + I2);

 return Tn;
}



