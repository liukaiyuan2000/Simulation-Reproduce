#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List DGP_cpp(int n, int Model_case, int H_case, int j_case){
  List res;
  if(Model_case == 1){
    NumericVector e(n), x1(n), x2(n), y(n);
    NumericMatrix X(n, 2);
    x1(0) = 0; 
    x2(0) = 0;
    for(int i=0; i<n; i++){
      e(i) = R::rnorm(0, 1);
      x1[i] = x1[i - 1]*0.5 + R::rnorm(0, 1);
      x2[i] = x2[i - 1]*0.5 + R::rnorm(0, 1);
      X(i, 0) = x1[i];
      X(i, 1) = x2[i];
    }
    double cn; 
    if(j_case == 1){
      cn = pow(n, -0.5)*sqrt(log(log(n)));
    }
    else{
      cn = pow(n, -7.0/18.0);
    }
    if(H_case == 1){//H0
      for(int i=0; i<n; i++){
        y(i) = x1(i) + x2(i) + e(i);
      }
    }
    else{//H1
      for(int i=0; i<n; i++){
        y(i) = x1(i) + x2(i) + e(i) + cn*(pow(x1(i), 2) + pow(x2(i), 2));
      }
    }
    res["X"] = X;
    res["y"] = y;
  }
  
  return res;
}

double kern(double t, int k) {
  // positive kernel functions
  // k is the switch between kernels
  double result;
  double pi = atan(1)*4;
  
  result = 0;
  if(k==1){ // Gaussian kernel
    result = sqrt(1/(2*pi))*exp(-pow(t,2)/2);
  }
  return result; 
}
double kern_conv1(double t, int k) {
  // positive kernel functions
  // k is the switch between kernels
  double result;
  double pi = atan(1)*4;
  
  result = 0;
  if(k==1){ // Gaussian kernel
    result = sqrt(1/(4*pi))*exp(-pow(t,2)/4);
  }
  return result; 
}
// [[Rcpp::export]]
double quantile_cpp(NumericVector x, double probs) {
  double quantiles;
  std::sort(x.begin(), x.end());
  
  double index = (x.size() - 1) * probs;
  int lowerIndex = floor(index);
  double weight = index - lowerIndex;
  quantiles = x[lowerIndex] * (1 - weight) + x[lowerIndex + 1] * weight;
  
  return quantiles;
}
double LL_x(NumericVector x, double h, int k) {
  int n = x.size();
  double LLi, res;
  res = 0.0;
  for (int i = 0; i < n; i ++){
    LLi = 0.0; 
    for (int j = 0; j < n; j ++){
      if(j != i){
        LLi = LLi + kern((x(i) - x(j))/h, k)/((n - 1)*h);
      }
    }
    //std::cout << LLi << std::endl;
    res = res + log(LLi)/n; 
  }
  
  return res;
}
double ISE_x(NumericVector x, double h, int k) {
  int n = x.size();
  double ISE1, ISE2, ISE;
  ISE1 = 0;
  ISE2 = 0;
  for (int i = 0; i < n; i ++){
    for (int j = 0; j < n; j ++){
      ISE1 += kern_conv1((x(i) - x(j))/h, k)/(pow(n, 2)*h);
      if(j != i){
        ISE2 += kern((x(i) - x(j))/h, k)/(n*(n - 1)*h);
      }
    }
  }
  ISE = ISE1 - 2*ISE2;
  return ISE;
}
NumericVector LCV_err(NumericVector x, NumericVector set_h, int k) {
  int num_bandwidths = set_h.size();
  NumericVector cv_error(num_bandwidths);
  
  for (int i = 0; i < num_bandwidths; i++) {
    double current_bandwidth = set_h[i];
    cv_error(i) = LL_x(x, current_bandwidth, k);
  }
  
  return cv_error;
}
NumericVector LSCV_err(NumericVector x, NumericVector set_h, int k) {
  int num_bandwidths = set_h.size();
  NumericVector cv_error(num_bandwidths);
  
  for (int i = 0; i < num_bandwidths; i++) {
    double current_bandwidth = set_h[i];
    cv_error(i) = ISE_x(x, current_bandwidth, k);
  }
  
  return cv_error;
}
// [[Rcpp::export]]
List ests_x_cpp(NumericMatrix X, NumericVector set_h, int k) {
  int n = X.rows();
  int d = X.cols();
  NumericMatrix kern_sum_mat(n, n);
  NumericVector bcv(d), pihat(n), pihat2(n);
  double v_2hat;
  
  for(int l=0; l<d; l++){
    NumericVector current_x(n);
    for(int i=0; i<n; i++){
      current_x(i) = X(i, l);
    }
    bcv(l) = set_h(which_min(LSCV_err(current_x, set_h, k)));
  }
  
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      kern_sum_mat(i, j) = 1.0;
      for(int l=0; l<d; l++){
        kern_sum_mat(i, j) = kern_sum_mat(i, j)*kern(abs(X(i,l) - X(j,l))/bcv(l), k);
      }
    }
  }
  double prod_bcv;
  prod_bcv = 1;
  for(int l=0; l<d; l++){
    prod_bcv *= bcv(l);
  }
  
  for(int i=0; i<n; i++){
    pihat(i) = 0.0; 
    for(int j=0; j<n; j++){
      pihat(i) += kern_sum_mat(i, j)/(prod_bcv*n);
    }
    pihat2(i) = pow(pihat(i), 2);
  }
  v_2hat = mean(pihat2);
  List res;
  res["kern_sum_mat"] = kern_sum_mat;
  res["pi.hat"] = pihat;
  res["v_2.hat"] = v_2hat;
  return res;
}
// [[Rcpp::export]]
List ests_yx_cpp(NumericVector y, NumericMatrix X) {
  int n = X.rows();
  int p = X.cols();
  Function solve("solve");
  
  NumericMatrix XtX(p, p);
  NumericVector Xty(p);
  for(int i = 0; i < n; i ++){
    for(int j = 0; j < p; j ++){
      for(int k = 0; k < p; k ++){
        XtX(j, k) += X(i, j) * X(i, k);
      }
      Xty(j) += X(i, j) * y(i);
    }
  }
  NumericVector thetahat = solve(XtX, Xty);
  NumericVector mhat(n), ehat(n);
  double mu_2hat;
  mu_2hat = 0;
  for(int i = 0; i < n; i ++){
    for(int j = 0; j < p; j ++){
      mhat(i) += X(i, j) * thetahat(j);
    }
    ehat(i) = y(i) - mhat(i);
    mu_2hat += pow(ehat(i), 2)/n;
  }
  List res;
  res["m.hat"] = mhat;
  res["e.hat"] = ehat;
  res["mu_2.hat"] = mu_2hat;
  return res;
}


// [[Rcpp::export]]
double h_ew_cpp(
    NumericVector y, NumericMatrix X, 
    List ests_x, List ests_yx, 
    double int_K_square, double Kern_conv3_zero
) {
  double d = X.ncol();
  int n = y.size();
  NumericVector pihat = as<NumericVector>(ests_x["pi.hat"]);
  NumericMatrix kern_sum_mat = as<NumericMatrix>(ests_x["kern_sum_mat"]);
  double v_2hat = as<double>(ests_x["v_2.hat"]);
  NumericVector ehat = as<NumericVector>(ests_yx["e.hat"]);
  double mu_2hat = as<double>(ests_yx["mu_2.hat"]);
  
  NumericVector delta_n(n), temp1(n), temp2(n), temp3(n);
  
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      temp1(i) = kern_sum_mat(i, j) + temp1(i);
      temp2(i) = kern_sum_mat(i, j)*ehat(j) + temp2(i);
    }
    delta_n(i) = temp2(i)/temp1(i);
    temp3(i) = pow(delta_n(i), 2)*pihat(i);
  }
  
  double C_n_square = Rcpp::mean(temp3)/(mu_2hat*sqrt(2*v_2hat*int_K_square));
  double t_nhat = n*C_n_square;
  double c_pihat = v_2hat/pow(Rcpp::mean(pihat), 3);
  double a_1hat = c_pihat*(pow(2, 0.5)*Kern_conv3_zero)/(3*pow(sqrt(int_K_square), 3));
  double h_ewhat = pow(a_1hat, (-1/(2*d)))*pow(t_nhat, (-3/(2*d)));
  return h_ewhat;
}

// [[Rcpp::export]]
double T_n_cpp(
    NumericVector y, NumericMatrix X, 
    List ests_x, List ests_yx, 
    double int_K_square, double Kern_conv3_zero
) {
  double d = X.ncol();
  int n = y.size();
  NumericVector pihat = as<NumericVector>(ests_x["pi.hat"]);
  NumericMatrix kern_sum_mat = as<NumericMatrix>(ests_x["kern_sum_mat"]);
  double v_2hat = as<double>(ests_x["v_2.hat"]);
  NumericVector ehat = as<NumericVector>(ests_yx["e.hat"]);
  double mu_2hat = as<double>(ests_yx["mu_2.hat"]);
  
  double sigma_n_square = 2*pow(mu_2hat, 2)*v_2hat*int_K_square;
  double temp1, temp2;
  temp1 = 0;
  temp2 = 0;
  for (int i = 0; i < n; i ++) {
    for (int j = 0; j < n; j ++) {
      temp1 = temp1 + ehat(i)*ehat(j)*kern_sum_mat(i, j);
    }
    temp1 = temp1 - ehat(i)*ehat(i)*kern_sum_mat(i, i);
  }
  
  temp2 = n*pow(h_ew_cpp(
    y, X, ests_x, ests_yx, int_K_square, Kern_conv3_zero
  ), d/2)*sqrt(sigma_n_square);
  
  //std::cout << temp2 << std::endl;
  double t_nhat = temp1/temp2;
  return t_nhat;
}

// [[Rcpp::export]]
double l_star_cpp(
    int B, double alpha, 
    NumericVector y, NumericMatrix X, 
    List ests_x, List ests_yx, 
    double int_K_square, double Kern_conv3_zero
) {
  int n = y.size();
  double mu_2hat = as<double>(ests_yx["mu_2.hat"]);
  NumericVector mhat = as<NumericVector>(ests_yx["m.hat"]);
  
  NumericVector T_n_star_hat_restore(B);
  for(int b = 0; b < B; b ++){
    NumericVector estar(n), ystar(n);
    for(int i = 0; i < n; i ++) {
      estar(i) = R::rnorm(0, 1);
      ystar(i) = mhat(i) + sqrt(mu_2hat)*estar(i);
    }
    T_n_star_hat_restore(b) = T_n_cpp(
      ystar, X, 
      ests_x, ests_yx_cpp(ystar, X), 
      int_K_square, Kern_conv3_zero
    );
  }
  
  double ln_star = quantile_cpp(T_n_star_hat_restore, alpha);
  return ln_star;
}


// [[Rcpp::export]]
double Sim_cpp(
    int n, int B, int nsim, int Model_case, int H_case, int j_case, int k, 
    double alpha, double int_K_square, double Kern_conv3_zero, 
    NumericVector set_h
) {
  NumericVector res(nsim);
  for(int i=0; i<nsim; i++){
    List data = DGP_cpp(n, Model_case, H_case, j_case);
    NumericVector y = as<NumericVector>(data["y"]);
    NumericMatrix X = as<NumericMatrix>(data["X"]);
    List ests_x_store = ests_x_cpp(X, set_h, k);
    List ests_yx_store = ests_yx_cpp(y, X);
    double current_l_star = l_star_cpp(B, alpha, y, X, ests_x_store, ests_yx_store, int_K_square, Kern_conv3_zero);
    double current_Tn = T_n_cpp(y, X, ests_x_store, ests_yx_store, int_K_square, Kern_conv3_zero);
    if(current_Tn > current_l_star){
      res(i) = 1.0;
    }
    else{
      res(i) = 0.0;
    }
    //std::cout << current_Tn << std::endl;
    //std::cout << current_l_star << std::endl;
    //std::cout << res << std::endl;
  }
  
  
  return mean(res);
}










