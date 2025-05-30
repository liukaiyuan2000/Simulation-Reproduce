#include <RcppEigen.h>
#include <cmath>

using namespace Rcpp;
using namespace Eigen;
using namespace std;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
MatrixXd aa(MatrixXd x, VectorXi idx) {
  return(x(idx, 0));
}
