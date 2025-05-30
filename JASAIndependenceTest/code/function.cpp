#include <RcppEigen.h>
#include <cmath>
#include <random>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;

double euclidean_distance(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2) {
  return std::sqrt((v1 - v2).squaredNorm());
}

double medmat(Eigen::MatrixXd mat) {
  NumericMatrix a = wrap(mat);
  return median(a);
}

Eigen::VectorXd sample(const Eigen::VectorXd& x, int size) {
  // Initialize random number generator
  std::random_device rd;
  std::mt19937 gen(rd());

  // Create a vector to store the sampled indices
  Eigen::VectorXi indices(size);

  // Generate random indices
  for (int i = 0; i < size; ++i) {
    std::uniform_int_distribution<> dis(0, x.size() - 1);
    indices(i) = dis(gen);
  }

  // Extract sampled elements
  Eigen::VectorXd sampled_elements(size);
  for (int i = 0; i < size; ++i) {
    sampled_elements(i) = x(indices(i));
  }

  return sampled_elements;
}

Eigen::VectorXd sample_cpp(const Eigen::VectorXd& x, int size) {
  return sample(x, size);
}

Eigen::VectorXi sample_n(int n, int k) {
  // Initialize a vector to store the sample
  Eigen::VectorXi sample(k);

  // Generate random indices without replacement
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 1);
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(indices.begin(), indices.end(), g);

  // Select the first k indices as the sample
  for (int i = 0; i < k; ++i) {
    sample[i] = indices[i];
  }

  return sample;
}

double sdmat(Eigen::MatrixXd mat) {
  NumericMatrix a = wrap(mat);
  return sd(a);
}

Eigen::MatrixXd computeRBFKernel(const Eigen::MatrixXd& X, double sigma) {
  const int n = X.rows();
  Eigen::MatrixXd K(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sq_dist = (X.row(i) - X.row(j)).squaredNorm();
      K(i, j) = std::exp(-sq_dist / (2 * pow(sigma, 2)));
      if (i != j)
        K(j, i) = K(i, j);
    }
  }
  return K;
}

Eigen::MatrixXd inchol_rcpp(const Eigen::MatrixXd& X, double sigma) {
  Eigen::MatrixXd K = computeRBFKernel(X, sigma);
  Eigen::LLT<Eigen::MatrixXd> lltOfK(K);
  Eigen::MatrixXd L = lltOfK.matrixL();
  return L;
}


// [[Rcpp::export]]
double slm_rcpp(Eigen::MatrixXd x, Eigen::MatrixXd y) {
  int n = x.rows();
  Eigen::MatrixXd x2vec(n*(n-1)/2, 1);
  int k = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      x2vec(k, 0) = pow(euclidean_distance(x.row(i), x.row(j)), 2);
      k = k + 1;
    }
  }
  double sigma = sqrt(0.5 * medmat(x2vec));
  // std::cout << sigma << std::endl;
  Eigen::MatrixXd K(n, n);
  Eigen::MatrixXd L(n, n);
  for (int i = 0; i < n; ++i) {
    K(i, i) = 1.0 / (std::sqrt(2.0 * M_PI)*sigma);
    L(i, i) = 1.0;
    for (int j = i + 1; j < n; ++j) {
      K(i, j) = (1.0 / (std::sqrt(2.0 * M_PI)*sigma)) * std::exp(-0.5 * std::pow(euclidean_distance(x.row(i), x.row(j))/sigma, 2));
      K(j, i) = K(i, j);
      L(i, j) = (euclidean_distance(y.row(i), y.row(j)) == 0) ? 1 : 0;
      L(j, i) = L(i, j);
    }
  }
  Eigen::VectorXd Lrowsum = L.rowwise().sum();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      L(i, j) = n * L(i, j) / Lrowsum(i) - 1;
    }
  }
  // std::cout << K << std::endl;
  // std::cout << L << std::endl;
  double res = (K*L).trace()/pow(n, 2)*sqrt(2*M_PI)*sigma;
  return res;
}


// [[Rcpp::export]]
Rcpp::List ptest_slm_rcpp(Eigen::MatrixXd x, Eigen::MatrixXd y, int B) {
  double Tn = slm_rcpp(x, y);

  // Simulate B samples and calculate p-value
  int count = 0;
  std::random_device rd;
  std::mt19937 gen(rd());
  for (int i = 0; i < B; ++i) {
    double h = slm_rcpp(x, sample_cpp(y, x.rows()));
    if (h >= Tn) {
      count++;
    }
  }
  double p = static_cast<double>(count + 1) / (B + 1);

  // Return results as a list
  return Rcpp::List::create(Rcpp::Named("ECCFIC") = Tn,
                            Rcpp::Named("pvalue") = p);
}


// [[Rcpp::export]]
double kem_rcpp(Eigen::MatrixXd x, Eigen::MatrixXd y) {

  int n = x.rows();
  int p = y.cols();

  Eigen::MatrixXd x2vec(n*(n-1)/2, 1);
  int k = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      x2vec(k, 0) = pow(euclidean_distance(x.row(i), x.row(j)), 2);
      k = k + 1;
    }
  }
  double sigma = sqrt(0.5 * medmat(x2vec));
  //std::cout << sigma << std::endl;
  double bw = 1.06 * sdmat(y) * pow(n, -1.0 / 5.0);
  Eigen::MatrixXd H = Eigen::MatrixXd::Identity(n, n) - Eigen::MatrixXd::Ones(n, n) / n;

  Eigen::MatrixXd tempx = inchol_rcpp(x, sigma);

  Eigen::JacobiSVD<Eigen::MatrixXd> svd_x(tempx, Eigen::ComputeThinU);
  Eigen::MatrixXd Ux = svd_x.matrixU();
  Eigen::VectorXd Sx2 = svd_x.singularValues().array().square();

  Eigen::MatrixXd tempy(n, n);
  for (int i = 0; i < n; ++i) {
    tempy(i, i) = 1.0 / (std::sqrt(2.0 * M_PI)*bw);
    for (int j = i + 1; j < n; ++j) {
      tempy(i, j) = (1.0 / (std::sqrt(2.0 * M_PI)*bw)) * std::exp(-0.5 * std::pow(euclidean_distance(y.row(i), y.row(j))/bw, 2));
      tempy(j, i) = tempy(i, j);
    }
  }
  tempy = tempy.array().colwise() / tempy.rowwise().sum().array();
  tempy = H * tempy.adjoint();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd_y(tempy, Eigen::ComputeThinU);
  Eigen::MatrixXd Uy = svd_y.matrixU();
  Eigen::VectorXd Sy2 = svd_y.singularValues().array().square();

  return (Ux * Sx2.asDiagonal() * Ux.transpose() * Uy * Sy2.asDiagonal() * Uy.transpose()).trace() / n;
}


// [[Rcpp::export]]
Rcpp::List ptest_kem_rcpp(Eigen::MatrixXd x, Eigen::MatrixXd y, int B) {
  int n = x.rows();
  double Tn = kem_rcpp(x, y);

  // Simulate B samples and calculate p-value
  int count = 0;
  std::random_device rd;
  std::mt19937 gen(rd());
  for (int i = 0; i < B; ++i) {
    Eigen::VectorXi idx = sample_n(n, n);
    double h = kem_rcpp(x, y(idx, 0));
    if (h >= Tn) {
      count++;
    }
  }
  double p = static_cast<double>(count + 1) / (B + 1);

  // Return results as a list
  return Rcpp::List::create(Rcpp::Named("ECCFIC") = Tn,
                            Rcpp::Named("pvalue") = p);
}

