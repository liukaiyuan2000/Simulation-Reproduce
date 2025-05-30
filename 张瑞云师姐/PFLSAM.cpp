#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

Eigen::MatrixXd matmatbind(Eigen::MatrixXd X, Eigen::MatrixXd Z) {
  int nrow = X.rows();
  int ncolX = X.cols();
  int ncolZ = Z.cols();
  Eigen::MatrixXd Y(nrow, ncolX + ncolZ);
  Y << X, Z;
  return Y;
}

Eigen::MatrixXd matvecbind(Eigen::MatrixXd X, Eigen::VectorXd Z) {
  int nrow = X.rows();
  int ncolX = X.cols();
  Eigen::MatrixXd Y(nrow, ncolX + 1);
  Y << X, Z;
  return Y;
}

Eigen::MatrixXd vecmatbind(Eigen::VectorXd X, Eigen::MatrixXd Z) {
  int nrow = Z.rows();
  int ncolZ = Z.cols();
  Eigen::MatrixXd Y(nrow, ncolZ + 1);
  Y << X, Z;
  return Y;
}

Eigen::MatrixXd vecvecbind(Eigen::VectorXd X, Eigen::VectorXd Z) {
  int nrow = X.rows();
  Eigen::MatrixXd Y(nrow, 2);
  Y << X, Z;
  return Y;
}

Eigen::MatrixXd ginvrcpp(Eigen::MatrixXd mat) {
  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(mat);
  return cod.pseudoInverse();
}

// [[Rcpp::export]]
List para_hat_rcpp(Eigen::MatrixXd Zn, Eigen::MatrixXd Wn, Eigen::VectorXd Yn, Eigen::MatrixXd Phi) {
  int n = Yn.size();

  // Calculate Qn
  Eigen::MatrixXd In = Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd WnYn = Wn * Yn;
  Eigen::MatrixXd Qn = matmatbind(WnYn, Zn);

  Eigen::MatrixXd P = Phi * ginvrcpp(Phi.transpose() * Phi) * Phi.transpose();
  Eigen::MatrixXd In_minus_P = In - P;

  // 2SLS
  // Step 1
  Eigen::MatrixXd Qnp = Qn * ginvrcpp(Qn.transpose() * Qn) * Qn.transpose();
  Eigen::MatrixXd temp1 = matmatbind(WnYn, matmatbind(Zn, Phi)).transpose();
  Eigen::MatrixXd temp2 = temp1 * Qnp;
  Eigen::VectorXd para_tilde = ginvrcpp(temp2 * temp1.transpose()) * temp2 * Yn;
  Eigen::MatrixXd H_tilde = matmatbind(Wn * ginvrcpp(In - para_tilde(0) * Wn) * matmatbind(Phi * para_tilde(3), Zn), Zn);

  // Step 2
  Eigen::MatrixXd M_tilde = H_tilde * ginvrcpp(H_tilde.transpose() * H_tilde) * H_tilde.transpose();
  Eigen::MatrixXd temp3 = Qn.transpose() * In_minus_P * M_tilde * In_minus_P;
  Eigen::MatrixXd temp4 = ginvrcpp(Phi.transpose() * Phi) * Phi.transpose();
  Eigen::VectorXd theta_bar = ginvrcpp(temp3 * Qn) * temp3 * Yn;
  Eigen::VectorXd alpha_bar = temp4 * (Yn - Qn * theta_bar);
  Eigen::MatrixXd H = matmatbind(Wn*ginvrcpp(In - theta_bar(0)*Wn)*(Phi * alpha_bar + Zn * theta_bar.tail(theta_bar.size() - 1)), Zn);
  Eigen::MatrixXd M = H * ginvrcpp(H.transpose() * H) * H.transpose();
  Eigen::MatrixXd temp5 = Qn.transpose() * In_minus_P * M * In_minus_P;
  Eigen::VectorXd theta_hat = ginvrcpp(temp5 * Qn) * temp5 * Yn;
  Eigen::VectorXd alpha_hat = temp4 * (Yn - Qn * theta_hat);

  return List::create(_["theta_hat"] = theta_hat, _["alpha_hat"] = alpha_hat);
}


// [[Rcpp::export]]
List CTn_cvm_rcpp(
    Eigen::MatrixXd Zn, Eigen::MatrixXd Wn, Eigen::VectorXd Yn, Eigen::MatrixXd Phi, Eigen::MatrixXd Xnt,
    Eigen::VectorXd theta_hat, Eigen::VectorXd alpha_hat, Eigen::VectorXd h_2_tilde
) {

  int time_num = Xnt.cols();
  double time_num_d = static_cast<double>(time_num);
  Eigen::VectorXd u = (Eigen::VectorXd::LinSpaced(time_num, 0.5, time_num - 0.5)) / time_num_d;
  int n = Yn.rows();
  Eigen::VectorXd varepsilon_hat = Yn - theta_hat[0] * Wn * Yn - Zn * theta_hat.tail(theta_hat.size() - 1) - Phi * alpha_hat;

  double h_mod = h_2_tilde.array().square().sum() + theta_hat.tail(theta_hat.size() - 1).array().square().sum();
  Eigen::VectorXd h = theta_hat.tail(theta_hat.size() - 1) / sqrt(h_mod);
  Eigen::VectorXd h_2 = h_2_tilde / sqrt(h_mod);

  Eigen::MatrixXd Uh = Zn * h + Xnt * h_2 / time_num_d;

  Eigen::VectorXd CRn_hat(time_num);
  for(int i = 0; i < time_num; i++) {
    double sum = 0;
    for(int j = 0; j < n; j++) {
      sum += varepsilon_hat[j] * (Uh(j) <= u(i) ? 1 : 0);
    }
    CRn_hat(i) = sum;
  }

  double CT_n = CRn_hat.array().square().mean();

  return List::create(_["varepsilon_hat"] = varepsilon_hat, _["CT_n"] = CT_n);
}

// [[Rcpp::export]]
List CTn_ks_rcpp(
    Eigen::MatrixXd Zn, Eigen::MatrixXd Wn, Eigen::VectorXd Yn, Eigen::MatrixXd Phi, Eigen::MatrixXd Xnt,
    Eigen::VectorXd theta_hat, Eigen::VectorXd alpha_hat, Eigen::VectorXd h_2_tilde
) {

  int time_num = Xnt.cols();
  double time_num_d = static_cast<double>(time_num);
  Eigen::VectorXd u = (Eigen::VectorXd::LinSpaced(time_num, 0.5, time_num - 0.5)) / time_num_d;
  int n = Yn.rows();
  Eigen::VectorXd varepsilon_hat = Yn - theta_hat[0] * Wn * Yn - Zn * theta_hat.tail(theta_hat.size() - 1) - Phi * alpha_hat;

  double h_mod = h_2_tilde.array().square().sum() + theta_hat.tail(theta_hat.size() - 1).array().square().sum();
  Eigen::VectorXd h = theta_hat.tail(theta_hat.size() - 1) / sqrt(h_mod);
  Eigen::VectorXd h_2 = h_2_tilde / sqrt(h_mod);

  Eigen::MatrixXd Uh = Zn * h + Xnt * h_2 / time_num_d;

  Eigen::VectorXd CRn_hat(time_num);
  for(int i = 0; i < time_num; i++) {
    double sum = 0;
    for(int j = 0; j < n; j++) {
      sum += varepsilon_hat[j] * (Uh(j) <= u(i) ? 1 : 0);
    }
    CRn_hat(i) = sum;
  }

  double CT_n = CRn_hat.maxCoeff();

  return List::create(_["varepsilon_hat"] = varepsilon_hat, _["CT_n"] = CT_n);
}


double bootstrap_single_cvm_rcpp(
    Eigen::MatrixXd Zn, Eigen::MatrixXd Wn, Eigen::MatrixXd Phi, Eigen::MatrixXd Xnt,Eigen::VectorXd theta_hat,
    Eigen::VectorXd alpha_hat, Eigen::VectorXd h_2_tilde, Eigen::VectorXd varepsilon_hat
) {

  int n = Xnt.rows();
  int time_num = Xnt.cols();
  Eigen::VectorXd u_star = (Eigen::VectorXd::LinSpaced(time_num, 0.5, time_num - 0.5)) / time_num;

  // Step 2
  Eigen::VectorXd varepsilon = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; i++) {
    varepsilon[i] = R::rnorm(0, 1);
  }
  Eigen::VectorXd varepsilon_star = varepsilon_hat.array() * varepsilon.array();

  Eigen::VectorXd Y_star = (Eigen::MatrixXd::Identity(n, n) - theta_hat[0] * Wn).inverse() * (Zn * theta_hat.tail(theta_hat.size() - 1) - Phi * alpha_hat + varepsilon_star);

  // Step 3
  List pavalue = para_hat_rcpp(Zn, Wn, Y_star, Phi);
  Eigen::VectorXd theta_hat_star = as<Eigen::VectorXd>(pavalue["theta_hat"]);
  Eigen::VectorXd alpha_hat_star = as<Eigen::VectorXd>(pavalue["alpha_hat"]);
  Eigen::VectorXd varepsilon_hat_star = Y_star - theta_hat_star[0] * Wn * Y_star - Zn * theta_hat_star.tail(theta_hat_star.size() - 1) - Phi * alpha_hat_star;

  // Step 4
  double h_mod = h_2_tilde.array().square().sum() + theta_hat.tail(theta_hat.size() - 1).array().square().sum();
  Eigen::VectorXd h = theta_hat.tail(theta_hat.size() - 1) / sqrt(h_mod);
  Eigen::VectorXd h_2 = h_2_tilde / sqrt(h_mod);

  Eigen::MatrixXd Uh = Zn * h + Xnt * h_2 / static_cast<double>(time_num);

  Eigen::VectorXd CRn_hat_star(time_num);
  for(int i = 0; i < time_num; i++) {
    double sum = 0;
    for(int j = 0; j < n; j++) {
      sum += varepsilon_hat_star[j] * (Uh(j) <= u_star(i) ? 1 : 0);
    }
    CRn_hat_star(i) = sum;
  }

  double CTn_star = CRn_hat_star.array().square().mean();

  return CTn_star;
}


double bootstrap_single_ks_rcpp(
    Eigen::MatrixXd Zn, Eigen::MatrixXd Wn, Eigen::MatrixXd Phi, Eigen::MatrixXd Xnt,Eigen::VectorXd theta_hat,
    Eigen::VectorXd alpha_hat, Eigen::VectorXd h_2_tilde, Eigen::VectorXd varepsilon_hat
) {

  int n = Xnt.rows();
  int time_num = Xnt.cols();
  Eigen::VectorXd u_star = (Eigen::VectorXd::LinSpaced(time_num, 0.5, time_num - 0.5)) / time_num;

  // Step 2
  Eigen::VectorXd varepsilon = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; i++) {
    varepsilon[i] = R::rnorm(0, 1);
  }
  Eigen::VectorXd varepsilon_star = varepsilon_hat.array() * varepsilon.array();

  Eigen::VectorXd Y_star = (Eigen::MatrixXd::Identity(n, n) - theta_hat[0] * Wn).inverse() * (Zn * theta_hat.tail(theta_hat.size() - 1) - Phi * alpha_hat + varepsilon_star);

  // Step 3
  List pavalue = para_hat_rcpp(Zn, Wn, Y_star, Phi);
  Eigen::VectorXd theta_hat_star = as<Eigen::VectorXd>(pavalue["theta_hat"]);
  Eigen::VectorXd alpha_hat_star = as<Eigen::VectorXd>(pavalue["alpha_hat"]);
  Eigen::VectorXd varepsilon_hat_star = Y_star - theta_hat_star[0] * Wn * Y_star - Zn * theta_hat_star.tail(theta_hat_star.size() - 1) - Phi * alpha_hat_star;

  // Step 4
  double h_mod = h_2_tilde.array().square().sum() + theta_hat.tail(theta_hat.size() - 1).array().square().sum();
  Eigen::VectorXd h = theta_hat.tail(theta_hat.size() - 1) / sqrt(h_mod);
  Eigen::VectorXd h_2 = h_2_tilde / sqrt(h_mod);

  Eigen::MatrixXd Uh = Zn * h + Xnt * h_2 / static_cast<double>(time_num);

  Eigen::VectorXd CRn_hat_star(time_num);
  for(int i = 0; i < time_num; i++) {
    double sum = 0;
    for(int j = 0; j < n; j++) {
      sum += varepsilon_hat_star[j] * (Uh(j) <= u_star(i) ? 1 : 0);
    }
    CRn_hat_star(i) = sum;
  }

  double CTn_star = CRn_hat_star.maxCoeff();

  return CTn_star;
}

// [[Rcpp::export]]
Eigen::VectorXd bootstrap_cvm_rcpp(
    Eigen::MatrixXd Zn, Eigen::MatrixXd Wn, Eigen::MatrixXd Phi, Eigen::MatrixXd Xnt,Eigen::VectorXd theta_hat,
    Eigen::VectorXd alpha_hat, Eigen::VectorXd h_2_tilde, Eigen::VectorXd varepsilon_hat, int B
) {
  Eigen::VectorXd CTn_star_sum = Eigen::VectorXd::Zero(B);
  for (int i = 0; i < B; i++) {
    CTn_star_sum[i] = bootstrap_single_cvm_rcpp(Zn, Wn, Phi, Xnt, theta_hat, alpha_hat, h_2_tilde, varepsilon_hat);
  }
  return CTn_star_sum;
}

// [[Rcpp::export]]
Eigen::VectorXd bootstrap_ks_rcpp(
    Eigen::MatrixXd Zn, Eigen::MatrixXd Wn, Eigen::MatrixXd Phi, Eigen::MatrixXd Xnt,Eigen::VectorXd theta_hat,
    Eigen::VectorXd alpha_hat, Eigen::VectorXd h_2_tilde, Eigen::VectorXd varepsilon_hat, int B
) {
  Eigen::VectorXd CTn_star_sum = Eigen::VectorXd::Zero(B);
  for (int i = 0; i < B; i++) {
    CTn_star_sum[i] = bootstrap_single_ks_rcpp(Zn, Wn, Phi, Xnt, theta_hat, alpha_hat, h_2_tilde, varepsilon_hat);
  }
  return CTn_star_sum;
}






