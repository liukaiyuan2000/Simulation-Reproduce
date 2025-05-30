#include <Rcpp.h>
using namespace Rcpp;

double kern(double t, int k) {
  // positive kernel functions
  // k is the switch between kernels
  double result;
  double pi = atan(1)*4;

  result = 0;
  if(k==1){ // uniform kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = 1;
  }
  if(k==2){ // triangular kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = (1-t)*2;
  }
  if(k==3){ // Epanechnikov kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = (1-pow(t,2))*3/2;
  }
  if(k==4){ // biweight kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = pow(1-pow(t,2),2)*15/8;
  }
  if(k==5){ // triweight kernel
    if ((t < 0)||(t > 1))
      result = 0;
    else
      result = pow(1-pow(t,2),3)*35/16;
  }
  if(k==6){ // Gaussian kernel
    if (t < 0)
      result = 0;
    else
      result = sqrt(2/pi)*exp(-pow(t,2)/2);
  }
  return result;
}

double rho_tau(double x, double tau){
  double res;
  if(x < 0){
    res = tau*x - x;
  }
  else{
    res = tau*x;
  }
  return res;
}

NumericVector rq_cpp(NumericVector Y, NumericMatrix X, double tau, NumericVector K) {
  // 创建一个数据框
  DataFrame data = DataFrame::create(Named("Y") = Y, Named("X") = X);

  // 设置tau和权重K的值
  NumericVector tau_vec = NumericVector::create(tau);
  NumericVector K_vec = K;

  // 加载quantreg包
  Environment quantreg = Environment("package:quantreg");
  Function rq = quantreg["rq"];

  // 创建formula字符串
  String formula = "Y ~ . - 1";

  // 调用rq函数
  List result = rq(Named("formula") = formula,
                   Named("data") = data,
                   Named("tau") = tau_vec,
                   Named("weights") = K_vec);
  NumericVector res = result["coefficients"];

  // 返回结果
  return res;
}


// Local Linear Quantile Regression for Rcpp
// [[Rcpp::export]]
NumericMatrix llqr_cpp(
    NumericMatrix iC, NumericMatrix iCnew, NumericMatrix D, NumericVector Y,
    NumericVector h, int kernI, double tau
) {
  int n = iC.nrow();
  int m = D.nrow();
  int J = iCnew.ncol();
  NumericMatrix C1(n, J);
  NumericVector K(n);
  NumericMatrix res(m, J);
  for (int k = 0; k < m; k++) {
    // for each function
    // create matrix with rows coefficients corresponding to X - x
    for (int i = 0; i < n; i++) {
      C1(i, 0) = 1;
      for (int j = 1; j < J; j++) {
        // C1[i, j] = C[i, j] - Cnew[k, j], i > 1
        C1(i, j) = iC(i, j) - iCnew(k, j);
      }
      K(i) = kern(D(k, i) / h(k), kernI);
    }
    NumericVector coef = rq_cpp(Y, C1, tau, K);
    for (int j = 0; j < J; j++) {
      res(k, j) = coef(j);
    }
  }

  return res;
}

// [[Rcpp::export]]
NumericMatrix llqr_cv_cpp(
    NumericMatrix C, NumericMatrix D, NumericVector Y,
    NumericMatrix H, int kernI, double tau
) {
  int i,j,k,l,hi;
  int n = D.nrow();
  int J = C.ncol();
  int nH = H.ncol();
  NumericMatrix Ci((n - 1), J), Di(1, (n - 1)), Cinew(1, J), Yt(1, J), CV(nH, J);
  NumericVector Yi(n - 1), h(1);

  // for each bandwidth h leave-one-out cross-validation
  for (hi = 0; hi < nH; hi++) {
    for(j = 0; j < J; j++){
      CV(hi, j) = 0;
    }
    for (i = 0; i < n; i++) {
      // cross-validation for first n functions only in a leave-one-out sense
      for (j = 0; j < J; j++) {
        Cinew(0, j) = C(i, j);
      }
      k = 0;
      for (j = 0; j < n; j++) {
        if (j != i) {
          Di(0, k) = D(i, j);
          Yi(k) = Y(j);
          for(l = 0; l < J; l++){
            Ci(k, l) = C(j, l);  // Ci[k,l] = C[j,l]
          }
          k++;
        }
      }
      h = H(i, hi);
      // Call llqr_cpp function
      Yt = llqr_cpp(Ci, Cinew, Di, Yi, h, kernI, tau);
      // cross-validate
      k = J - 1;
      double temp = Y(i)-Yt(0,0);
      CV(hi, k) += pow(temp, 2);
    }
    j = J - 1;
    CV(hi, j) = CV(hi, j) / n;
  }
  return CV;
}


// [[Rcpp::export]]
NumericMatrix llqr_leave_cpp(
    NumericMatrix iC, NumericMatrix iCnew, NumericMatrix D, NumericVector Y,
    NumericVector h, int kernI, double tau
) {
  int n = iC.nrow();
  int m = D.nrow();
  int J = iCnew.ncol();
  NumericMatrix C1(n, J);
  NumericVector K(n);
  NumericMatrix res(m, J);
  for (int k = 0; k < m; k++) {
    // for each function
    // create matrix with rows coefficients corresponding to X - x
    for (int i = 0; i < n; i++) {
      C1(i, 0) = 1;
      for (int j = 1; j < J; j++) {
        // C1[i, j] = C[i, j] - Cnew[k, j], i > 1
        C1(i, j) = iC(i, j) - iCnew(k, j);
      }
      if(i == k){
        K(i) = 0;
      }
      else{
        K(i) = kern(D(k, i) / h(k), kernI);
      }
    }
    NumericVector coef = rq_cpp(Y, C1, tau, K);
    for (int j = 0; j < J; j++) {
      res(k, j) = coef(j);
    }
  }
  return res;
}













