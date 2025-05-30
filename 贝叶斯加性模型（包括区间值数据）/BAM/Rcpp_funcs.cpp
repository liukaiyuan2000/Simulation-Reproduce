#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <random>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace Eigen;
using namespace std;

NumericMatrix mvrnorm_rcpp(int n, NumericVector mu, NumericMatrix Sigma) {
  int p = mu.size();
  std::mt19937 generator(std::random_device{}());
  std::normal_distribution<double> distribution(0.0, 1.0);
  
  NumericMatrix result(n, p);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      double sum = 0.0;
      for (int k = 0; k < p; k++) {
        sum += Sigma(j, k) * distribution(generator);
      }
      result(i, j) = mu[j] + sum;
    }
  }
  
  return result;
}

NumericVector rgamma_rcpp(int n, double shape, double scale) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::gamma_distribution<double> gamma(shape, scale);
  
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = gamma(gen);
  }
  
  return result;
}

// [[Rcpp::export]]
double matdet(Eigen::MatrixXd mat) {
  return mat.determinant();
}


Eigen::MatrixXd eigenmattrans(NumericMatrix X) {
  
  Eigen::Map<Eigen::MatrixXd> A(as<Eigen::Map<Eigen::MatrixXd> >(X));
  
  return(A);
}

Eigen::MatrixXd vec_to_rowmat(NumericVector vec) {
  Eigen::Map<Eigen::MatrixXd> mat(vec.begin(), 1, vec.size());
  
  return mat;
}

NumericMatrix rcppmattrans(Eigen::MatrixXd X) {
  
  SEXP s = wrap(X);
  NumericMatrix w(s);
  
  return(w);
}

Eigen::MatrixXd matinv(NumericMatrix X) {
  Eigen::MatrixXd mat = eigenmattrans(X);
  return mat.inverse();
}

NumericMatrix center_rcpp(NumericMatrix x) {
  int n = x.nrow(), p = x.ncol();
  
  NumericVector colMeans(p);
  
  for (int j = 0; j < p; j++) {
    NumericVector column = x(_, j);
    colMeans[j] = mean(column);
  }
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      x(i, j) = (x(i, j) - colMeans[j]);
    }
  }
  
  return x;
}


NumericMatrix mattrans(NumericMatrix mat) {
  int rows = mat.nrow();
  int cols = mat.ncol();
  
  NumericMatrix result(cols, rows);
  
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      result(j, i) = mat(i, j);
    }
  }
  
  return result;
}


NumericMatrix matmatmult(NumericMatrix A, NumericMatrix B) {
  Eigen::MatrixXd mA = eigenmattrans(A);
  Eigen::MatrixXd mB = eigenmattrans(B);
  
  return Rcpp::wrap(mA*mB);
}

NumericVector matvecmult(NumericMatrix A, NumericVector v) {
  int nrow_A = A.nrow();
  int ncol_A = A.ncol();
  
  NumericVector result(nrow_A);
  
  for (int i = 0; i < nrow_A; i++) {
    for (int j = 0; j < ncol_A; j++) {
      result[i] += A(i, j) * v[j];
    }
  }
  
  return result;
}

NumericMatrix matmatadd(NumericMatrix A, NumericMatrix B) {
  if (A.nrow() != B.nrow() || A.ncol() != B.ncol()) {
    throw std::range_error("Matrices must have the same dimensions");
  }
  
  NumericMatrix C(A.nrow(), A.ncol());
  
  for (int i = 0; i < A.nrow(); i++) {
    for (int j = 0; j < A.ncol(); j++) {
      C(i, j) = A(i, j) + B(i, j);
    }
  }
  
  return C;
}

// Function to generate runif
NumericVector runif_rcpp(int n, double min, double max) {
  NumericVector res(n);
  for (int i = 0; i < n; i++) {
    res[i] = R::runif(min, max);
  }
  return res;
}
// [[Rcpp::export]]
List DGP_rcpp(int n) {
  Rcpp::Function bs("bs");
  NumericVector x1 = runif_rcpp(n, -1, 1);
  NumericVector x2 = runif_rcpp(n, 0, 1);
  NumericVector x3 = runif_rcpp(n, 0, 1);
  NumericMatrix x(n, 3);
  x(_, 0) = x1;
  x(_, 1) = x2;
  x(_, 2) = x3;
  NumericMatrix z(n, 2);
  for(int i = 0; i < n; i++){
    z(i, 0) = 1.0;
  }
  z(_, 1) = runif_rcpp(n, -1, 1);
  
  NumericVector gamma = NumericVector::create(-0.5, 0.5); // Assume gamma is a vector of length 2
  double mu = 1;
  NumericVector sig2 = exp(matvecmult(z, gamma));
  NumericVector eps = rnorm(n) * sqrt(sig2);
  NumericVector y = mu + x1 + sin(2 * M_PI * x2) + cos(2 * M_PI * x3) + eps;
  int p = x.ncol();
  
  List B(p);
  for (int i = 0; i < p; i++) {
    NumericVector xi = x(_, i);
    NumericVector knots = NumericVector::create(mean(xi));
    NumericVector bs_func = bs(Rcpp::_["x"] = xi, Rcpp::_["degree"] = 3, Rcpp::_["knots"] = knots);
    B[i] = center_rcpp(as<NumericMatrix>(bs_func));
  }
  
  return List::create(
    Named("x") = x,
    Named("y") = y,
    Named("z") = z,
    Named("B") = B
  );
}

// [[Rcpp::export]]
double p_gamma_rcpp(
    NumericVector gamma, List alpha, double mu, NumericMatrix x,
    NumericVector y, NumericMatrix z, List B,
    NumericMatrix B_gamma, NumericVector gamma0, int K
) {
  int n = x.nrow();
  int p = B.size();
  NumericMatrix B_gamma_inv = rcppmattrans(matinv(B_gamma));
  NumericMatrix Sig(n, n);
  for (int i = 0; i < n; i++) {
    Sig(i, i) = exp(sum(z(i, _) * gamma));
  }
  NumericMatrix Sig_inv = rcppmattrans(matinv(Sig));
  NumericVector temp1(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++){
      NumericMatrix Bi = B[j];
      NumericVector alphai = alpha[j];
      for (int k = 0; k < K; k++) {
        temp1(i) += Bi(i, k) * alphai(k);
      }
    }
  }
  
  double f1;
  f1 = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++){
      f1 += (y(i) - temp1(i) - mu)*(y(j) - temp1(j) - mu)*Sig_inv(i, j);
    }
  }
  double f2;
  f2 = 0;
  for (int i = 0; i < p - 1; i++) {
    for (int j = 0; j < p - 1; j++){
      f2 += (gamma(i) - gamma0(i))*(gamma(j) - gamma0(j))*B_gamma_inv(i, j);
    }
  }
  double f = exp(-0.5 * (f1 + f2));
  double det_Sig = matdet(eigenmattrans(Sig));
  return f / pow(det_Sig, 0.5);
}
// [[Rcpp::export]]
List Gibbs_sample_rcpp(int n, int burn_in, int mc_num) {
  Rcpp::Function rnorm("rnorm");
  Rcpp::Function diag("diag");
  
  // Define prior parameters
  NumericVector gamma0 = NumericVector::create(-0.5, 0.5);
  double mu0 = 0;
  NumericMatrix B_gamma = diag(2);
  NumericMatrix B_gamma_inv = rcppmattrans(matinv(B_gamma));
  double a_tau = 1;
  double b_tau = 1;
  int K = 4;
  
  // Generate data using DGP function
  List data = DGP_rcpp(n);
  NumericMatrix x = as<NumericMatrix>(data["x"]);
  NumericMatrix z = as<NumericMatrix>(data["z"]);
  NumericVector y = as<NumericVector>(data["y"]);
  List B = data["B"];
  
  int p = B.size();
  int N = burn_in + mc_num;
  double mu_old = mu0;
  NumericVector gamma_old = gamma0;
  NumericVector gamma_new(2, 0.0);
  NumericMatrix Sigma_old = diag(exp(matvecmult(z, gamma_old)));
  
  NumericMatrix gamma_process(N, 2);
  NumericVector mu_process(N);
  List alpha_process(p);
  List alpha_iter(p);
  List alpha0(p);
  NumericVector p_acc(N);
  
  for (int i = 0; i < p; i++) {
    alpha_process[i] = NumericMatrix(N, K);
    alpha0[i] = NumericVector(K, 0.0);
  }
  alpha_iter = alpha0;
  
  for (int s = 0; s < N; s++) {
    Rcpp::Rcout << s << "\r";
    NumericMatrix Sigma_old_inv = rcppmattrans(matinv(Sigma_old));
    NumericVector z_gamma_old_mult = matvecmult(z, gamma_old);
    
    // calculate tau2.new
    List tau2_new(p);
    for (int i = 0; i < p; i++) {
      NumericVector alpha_iter_i = alpha_iter[i];
      NumericVector alpha0_i = alpha0[i];
      double tau2_temp = 0;
      for(int l = 0; l < K; l++){
        tau2_temp += pow(alpha_iter_i[l] - alpha0_i[l], 2);
      }
      tau2_new[i] = 1.0 / as<double>(rgamma_rcpp(1, K/2 + a_tau, 0.5 * (tau2_temp + 2 * b_tau)));
    }
    
    // calculate alpha.new
    for (int i = 0; i < p; i++) {
      NumericVector temp1(n, 0.0);
      for (int j = 0; j < p; j++){
        if(j != i){
          NumericMatrix B_j = B[j];
          NumericVector alpha_iter_j = alpha_iter[j];
          temp1 += matvecmult(B_j, alpha_iter_j);
        }
      }
      NumericMatrix b_alpha_star_inv(K, K);
      NumericMatrix t1 = diag(K);
      NumericMatrix t2 = matmatmult(matmatmult(mattrans(as<NumericMatrix>(B[i])), Sigma_old_inv), B[i]);
      for(int l = 0; l < K; l++){
        double tau2_new_temp = tau2_new[i];
        t1(l, l) = 1.0 / tau2_new_temp;
        b_alpha_star_inv(_, l) = t1(_, l) + t2(_, l);
      }
      NumericMatrix b_alpha_star = rcppmattrans(matinv(b_alpha_star_inv));
      NumericVector alpha0_star = matvecmult(
        b_alpha_star, matvecmult(t1, alpha0[i]) + matvecmult(matmatmult(mattrans(as<NumericMatrix>(B[i])), Sigma_old_inv), (y - mu_old - temp1))
      );
      alpha_iter[i] = as<NumericVector>(mvrnorm_rcpp(1, alpha0_star, b_alpha_star));
      NumericMatrix alpha_store_mat = alpha_process[i];
      alpha_store_mat(s,_) = as<NumericVector>(alpha_iter[i]);
    }
    
    // calculate mu.new
    NumericVector temp2(n, 0.0);
    for (int j = 0; j < p; j++) {
      temp2 += matvecmult(B[j], alpha_iter[j]);
    }
    
    NumericMatrix Sigma_old_inv_diag_mat = rcppmattrans(vec_to_rowmat(as<NumericVector>(diag(Sigma_old_inv))));
    double mu_sig = 1 / sum(Sigma_old_inv_diag_mat);
    NumericVector mu_mu = matvecmult(Sigma_old_inv_diag_mat, (y - temp2)) * mu_sig;
    double mu_new = as<double>(rnorm(1, mu_mu, sqrt(mu_sig)));
    
    // calculate gamma.new, MH Algorithm
    NumericMatrix Omega_gamma = 0.5 * matmatadd(matmatmult(
      matmatmult(mattrans(z), as<NumericMatrix>(diag(pow((y - mu_new - temp2), 2) / exp(z_gamma_old_mult)))), z
    ), B_gamma_inv);
    double sig_gamma2 = 1.5; // Example value, replace with actual value
    NumericVector mc_new = as<NumericVector>(mvrnorm_rcpp(1, gamma_old, sig_gamma2 * rcppmattrans(matinv(Omega_gamma))));
    p_acc[s] = std::min(
      1.0,
      p_gamma_rcpp(mc_new, alpha_iter, mu_new, x, y, z, B, B_gamma, gamma0, K) /
        p_gamma_rcpp(gamma_old, alpha_iter, mu_new, x, y, z, B, B_gamma, gamma0, K)
    );
    if(as<double>(runif_rcpp(1, 0, 1)) < p_acc[s]){
      gamma_new = mc_new;
    }
    else{
      gamma_new = gamma_old;
    }
    gamma_process(s, _) = gamma_old;
    mu_process[s] = mu_old;
    
    gamma_old = gamma_new;
    Sigma_old = diag(exp(matvecmult(z, gamma_new)));
    mu_old = mu_new;
    
    
  }
  
  return List::create(
    Named("data") = data,
    Named("gamma.process") = gamma_process,
    Named("alpha.process") = alpha_process,
    Named("mu.process") = mu_process
  );
}


