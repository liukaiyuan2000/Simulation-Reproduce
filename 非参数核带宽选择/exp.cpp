#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double quantileExample(NumericVector x, double probs) {
  double quantiles;
  std::sort(x.begin(), x.end());
  
  double index = (x.size() - 1) * probs;
  int lowerIndex = floor(index);
  double weight = index - lowerIndex;
  quantiles = x[lowerIndex] * (1 - weight) + x[lowerIndex + 1] * weight;
  
  return quantiles;
}
