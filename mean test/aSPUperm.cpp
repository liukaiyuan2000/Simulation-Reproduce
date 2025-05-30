#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <random>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// Function to compute covariance matrix
Eigen::MatrixXd cov_rcpp(const Eigen::MatrixXd& X) {
  Eigen::MatrixXd centered = X.rowwise() - X.colwise().mean();
  return (centered.adjoint() * centered) / double(X.rows() - 1);
}

// Function to compute ranks
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

// Function to compute abs of a NumericVector
Rcpp::NumericVector abs_vector(Rcpp::NumericVector x) {
  Rcpp::NumericVector result(x.size());
  std::transform(x.begin(), x.end(), result.begin(), [](double val) { return std::abs(val); });
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector aSPUperm_rcpp(
    const Eigen::MatrixXd& sam1, const Eigen::MatrixXd& sam2, 
    Rcpp::NumericVector pow = Rcpp::NumericVector::create(1, 2, 3, 4, 5, 6, 10), 
    int n_perm = 1000
) {
  int n1 = sam1.rows();
  int n2 = sam2.rows();
  int p = sam1.cols();
  
  Eigen::MatrixXd Sn = ((n1 - 1) * cov_rcpp(sam1) + (n2 - 1) * cov_rcpp(sam2)) / (n1 + n2 - 2);
  Eigen::MatrixXd sam(n1 + n2, p);
  sam << sam1, sam2;
  
  Eigen::VectorXd diff = sam1.colwise().mean() - sam2.colwise().mean();
  Rcpp::NumericVector Ts(pow.size());
  
  for (int j = 0; j < pow.size(); ++j) {
    if (pow[j] < 10) {
      Ts[j] = (diff.array().pow(pow[j])).sum();
    } else {
      Eigen::VectorXd diag_Sn = Sn.diagonal();
      for (int i = 0; i < diag_Sn.size(); ++i) {
        if (diag_Sn[i] <= 1e-10) diag_Sn[i] = 1e-10;
      }
      Ts[j] = (diff.array().square() / diag_Sn.array()).maxCoeff();
    }
  }
  // std::cout << Ts << std::endl;
  Rcpp::NumericVector p_spu(pow.size(), 0.0);
  Eigen::MatrixXd Ts_perm(pow.size(), n_perm);
  
  std::random_device rd;
  std::mt19937 g(rd());
  
  for (int b = 0; b < n_perm; ++b) {
    Eigen::VectorXi perm = Eigen::VectorXi::LinSpaced(n1 + n2, 0, n1 + n2 - 1);
    std::shuffle(perm.data(), perm.data() + perm.size(), g);
    
    Eigen::MatrixXd sam_perm(n1 + n2, p);
    for (int i = 0; i < n1 + n2; ++i) {
      sam_perm.row(i) = sam.row(perm[i]);
    }
    
    Eigen::MatrixXd sam1_perm = sam_perm.topRows(n1);
    Eigen::MatrixXd sam2_perm = sam_perm.bottomRows(n2);
    
    Eigen::MatrixXd Sn_perm = ((n1 - 1) * cov_rcpp(sam1_perm) + (n2 - 1) * cov_rcpp(sam2_perm)) / (n1 + n2 - 2);
    Eigen::VectorXd diff_perm = sam1_perm.colwise().mean() - sam2_perm.colwise().mean();
    
    for (int j = 0; j < pow.size(); ++j) {
      if (pow[j] < 10) {
        Ts_perm(j, b) = (diff_perm.array().pow(pow[j])).sum();
      }
      if (pow[j] == 10) {
        Eigen::VectorXd diag_Sn_perm = Sn_perm.diagonal();
        for (int i = 0; i < diag_Sn_perm.size(); ++i) {
          if (diag_Sn_perm[i] <= 1e-10) diag_Sn_perm[i] = 1e-10;
        }
        Ts_perm(j, b) = (diff_perm.array().square() / diag_Sn_perm.array()).maxCoeff();
      }
    }
  }
  
  Rcpp::NumericVector minp_perm(n_perm);

  for (int j = 0; j < pow.size(); ++j) {
    p_spu[j] = (
      std::count_if(
        Ts_perm.row(j).begin(), Ts_perm.row(j).end(),
        [Ts, j](double val) { return std::abs(val) >= std::abs(Ts[j]); }) + 1.0
    ) / (n_perm + 1.0);

    Rcpp::NumericVector p_spu_perm = (n_perm + 1.0 - rank_rcpp(abs_vector(Rcpp::wrap(Ts_perm.row(j))))) / double(n_perm);
    if (j == 0) {
      minp_perm = p_spu_perm;
    } else {
      for (int i = 0; i < n_perm; ++i) {
        if (minp_perm[i] > p_spu_perm[i]) {
          minp_perm[i] = p_spu_perm[i];
        }
      }
    }
  }
  
  double min_p_spu = Rcpp::min(p_spu);
  double p_aspu = (
    std::count_if(minp_perm.begin(), minp_perm.end(), 
                  [min_p_spu](double val) { return val <= min_p_spu; }) + 1
  ) / double(n_perm + 1);
  
  Rcpp::NumericVector pvs(pow.size() + 1);
  std::copy(p_spu.begin(), p_spu.end(), pvs.begin());
  pvs[pow.size()] = p_aspu;
  
  // Rcpp::CharacterVector names(pow.size() + 1);
  // for (int i = 0; i < pow.size() - 1; ++i) {
  //   names[i] = "SPU_" + std::to_string(int(pow[i]));
  // }
  // names[pow.size() - 1] = "SPU_Inf";
  // names[pow.size()] = "aSPU";
  // pvs.attr("names") = names;
  return pvs;
}













