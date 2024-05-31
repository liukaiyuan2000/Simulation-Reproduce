#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix M_boot_func(int Nk_boot, NumericVector X, NumericVector Y, int n, double epsilon, NumericVector par_c) {
  Rcpp::Environment my_Env = Rcpp::Environment::global_env();
  Rcpp::Function M_tests_boot_R = my_Env["M.tests.boot"];
  NumericMatrix val_M_boot(2, Nk_boot);
  
  for (int j = 0; j < Nk_boot; j++) {
    
    NumericVector Xi = rnorm(n);
    
    NumericVector outcome = M_tests_boot_R(X, Y, n, epsilon, Xi, par_c);
    
    val_M_boot(0, j) = outcome[0];
    val_M_boot(1, j) = outcome[3];
  }
  
  return val_M_boot;
}



