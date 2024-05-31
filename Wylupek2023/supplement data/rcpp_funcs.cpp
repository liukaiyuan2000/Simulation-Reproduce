#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ecdf_rcpp(NumericVector x, NumericVector y) {
  int n = y.size();
  int m = x.size();
  NumericVector result(m);
  
  for (int i = 0; i < m; i++) {
    result[i] = sum(y <= x[i]) / double(n);
  }
  
  return result;
}
// [[Rcpp::export]]
NumericVector sort_c(NumericVector x){
  std::sort(x.begin(), x.end());
  return(x);
}

// [[Rcpp::export]]
NumericVector M_tests(NumericVector X, NumericVector Y, int n, double epsilon = 0.0001) {
  int m = X.size();
  int p = Y.size();
  NumericVector Z;
  
  for (int i = 0; i < m; i++) {
    Z.push_back(X[i]);
  }
  for (int i = 0; i < p; i++) {
    Z.push_back(Y[i]);
  }
  
  NumericVector F = ecdf_rcpp(Z, X); 
  NumericVector G = ecdf_rcpp(Z, Y); 
  NumericVector H = ecdf_rcpp(Z, pmax(X, Y)); 
  NumericVector E = F - G;
  double KS = sqrt(n) * max(abs(E));
  
  NumericVector w = sqrt(abs(F + G - 2 * H - pow((F - G), 2)));
  double wKS = sqrt(n) * max(abs(ifelse(w == 0, 0, E / w)));
  
  NumericVector Z_t;
  NumericVector sortZ = sort_c(Z);
  for (int i = std::ceil(epsilon * 2 * n) - 1; i < std::floor((1 - epsilon) * 2 * n) - 1; i++) {
    Z_t.push_back(sortZ[i]);
  }
  NumericVector F_t = ecdf_rcpp(Z_t, X); 
  NumericVector G_t = ecdf_rcpp(Z_t, Y); 
  NumericVector H_t = ecdf_rcpp(Z_t, pmax(X, Y)); 
  NumericVector E_t = F_t - G_t;
  NumericVector w_t = sqrt(abs(F_t + G_t - 2 * H_t - pow(F_t - G_t, 2)));
  double WKS = sqrt(n) * max(abs(ifelse(w_t == 0, 0, E_t / w_t)));
  
  double M = std::max(KS, WKS);
  
  NumericVector all = NumericVector::create(KS, wKS, WKS, M);
  return all;
}


// [[Rcpp::export]]
NumericVector M_tests_boot(
    NumericVector X, NumericVector Y, int n, 
    NumericVector vec_Xi, double par_c, double epsilon=0.0001
) {
  
  
  NumericVector Z;
  
  for (int i = 0; i < X.size(); i++) {
    Z.push_back(X[i]);
  }
  for (int i = 0; i < Y.size(); i++) {
    Z.push_back(Y[i]);
  }
  NumericVector Z_s = sort_c(Z);
  
  int m = 2*n;
  
  NumericVector E_boot(m);
  NumericVector w_boot(m);
  
  for (int i = 0; i < m; i++) {
    double sum_X = 0; 
    double sum_w1 = 0; 
    double sum_w2 = 0; 
    for(int j = 0; j < n; j++) {
      sum_X += vec_Xi[j] * (((X[j] <= Z_s[i]) ? 1 : 0) - ((Y[j] <= Z_s[i]) ? 1 : 0)) / sqrt(n);
      sum_w1 += (abs(vec_Xi[j] / sqrt(2 / M_PI)) * (((X[j] <= Z_s[i]) ? 1 : 0) + 
        ((Y[j] <= Z_s[i]) ? 1 : 0) - 2 * ((std::max(X[j], Y[j]) <= Z_s[i]) ? 1 : 0))) / n;
      sum_w2 += (abs(vec_Xi[j] / sqrt(2 / M_PI)) * 
        (((X[j] <= Z_s[i]) ? 1 : 0) - ((Y[j] <= Z_s[i]) ? 1 : 0)));
    }
    E_boot[i] = sum_X;
    w_boot[i] = sum_w1 - pow(sum_w2, 2) / pow(n, 2);
  }
  
  double KS_boot = max(abs(E_boot));
  double wKS_boot = max(abs(ifelse(w_boot == 0, 0, E_boot / sqrt(abs(par_c * w_boot)))));
  
  
  NumericVector Z_t;
  for (int i = std::ceil(epsilon * m) - 1; i < std::floor((1 - epsilon) * m) - 1; i++) {
    Z_t.push_back(Z_s[i]);
  }
  int dl = Z_t.size();
  
  NumericVector E_t_boot(dl);
  NumericVector w_t_boot(dl);
  
  for (int i = 0; i < dl; i++) {
    double sum_X_t = 0; 
    double sum_w_t1 = 0; 
    double sum_w_t2 = 0; 
    for(int j = 0; j < n; j++) {
      sum_X_t += vec_Xi[j] * (((X[j] <= Z_t[i]) ? 1 : 0) - ((Y[j] <= Z_t[i]) ? 1 : 0)) / sqrt(n);
      sum_w_t1 += (abs(vec_Xi[j] / sqrt(2 / M_PI)) * (((X[j] <= Z_t[i]) ? 1 : 0) + 
        ((Y[j] <= Z_t[i]) ? 1 : 0) - 2 * ((std::max(X[j], Y[j]) <= Z_t[i]) ? 1 : 0))) / n;
      sum_w_t2 += (abs(vec_Xi[j] / sqrt(2 / M_PI)) * 
        (((X[j] <= Z_t[i]) ? 1 : 0) - ((Y[j] <= Z_t[i]) ? 1 : 0)));
    }
    E_t_boot[i] = sum_X_t;
    w_t_boot[i] = sum_w_t1 - pow(sum_w_t2, 2) / pow(n, 2);
  }
  
  double WKS_boot = max(abs(ifelse(w_t_boot == 0, 0, E_t_boot / sqrt(abs(par_c * w_t_boot)))));
  double M_boot = std::max(KS_boot, WKS_boot);
  
  NumericVector all_boot = NumericVector::create(KS_boot, wKS_boot, WKS_boot, M_boot);
  return all_boot;
}


// [[Rcpp::export]]
NumericMatrix M_boot_func(int Nk_boot, NumericVector X, NumericVector Y, int n, double epsilon, double par_c) {
  NumericMatrix val_M_boot(2, Nk_boot);
  
  for (int j = 0; j < Nk_boot; j++) {
    
    NumericVector Xi = rnorm(n);
    
    NumericVector outcome = M_tests_boot(X, Y, n, Xi, par_c, epsilon);
    
    val_M_boot(0, j) = outcome[0];
    val_M_boot(1, j) = outcome[3];
  }
  
  return val_M_boot;
}


