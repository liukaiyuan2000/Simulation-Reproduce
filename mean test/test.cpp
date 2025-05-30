#include <RcppEigen.h>
#include <algorithm> // for std::transform
#include <cmath>     // for std::abs
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Rcpp::NumericVector rank_rcpp(Rcpp::NumericVector x) {
  int n = x.size();
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  
  std::sort(indices.begin(), indices.end(), [&x](int i1, int i2) { return x[i1] < x[i2]; });
  
  Rcpp::NumericVector ranks(n);
  for (int i = 0; i < n; ++i) {
    ranks[indices[i]] = i + 1; // ranks start from 1
  }
  
  return ranks;
}

// [[Rcpp::export]]
double aa(NumericVector x) {
  
  return min(x);
}


