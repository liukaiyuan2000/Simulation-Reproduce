#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double fkernel(double d, double h) {
  return std::max(0.0, 0.75 * (1 - std::pow(d / h, 2)));
}
// Function to calculate the inner product with delta scaling
// [[Rcpp::export]]
double inpro(const VectorXd &x, const VectorXd &y, const VectorXd &t) {
  int n = x.size();
  double delta = (t(n - 1) - t(0)) / (n - 1);
  return x.dot(y) * delta;
}

// [[Rcpp::export]]
List pca_rcppeigen(const Eigen::MatrixXd &X, int m, const Eigen::VectorXd &t) {
  int n = X.rows();
  
  // Eigen decomposition
  Eigen::MatrixXd cov_mat = (X.transpose() * X) / (n - 1); // calculate covariance
  SelfAdjointEigenSolver<MatrixXd> eig_solver(cov_mat);
  MatrixXd eig_vectors = eig_solver.eigenvectors().rightCols(m);
  
  // Scaling eigenvectors by sqrt(t.size())
  double scale_factor = sqrt(t.size());
  eig_vectors *= scale_factor;
  
  // Eigen scores and values initialization
  MatrixXd eig_scores(n, m);
  VectorXd eig_values(m);
  
  if (m > 1) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        eig_scores(i, j) = inpro(X.row(i), eig_vectors.col(j), t);
      }
    }
    eig_values = eig_scores.array().square().colwise().mean();
  } else {
    VectorXd eig_scores_single(n);
    for (int i = 0; i < n; i++) {
      eig_scores_single(i) = inpro(X.row(i), eig_vectors.col(0), t);
    }
    eig_values(0) = eig_scores_single.array().square().mean();
    eig_scores.col(0) = eig_scores_single;
  }
  
  return List::create(
    Named("eig.values") = eig_values,
    Named("eig.functions") = eig_vectors,
    Named("h") = eig_scores
  );
}

// [[Rcpp::export]]
List BE_stat(const Eigen::MatrixXd &X, const Eigen::VectorXd &Y, const Eigen::MatrixXd &Z, int m, const Eigen::VectorXd &tm) {
  int n = X.rows();
  
  // PCA results
  List pca_res = pca_rcppeigen(X, m, tm);
  MatrixXd H = as<MatrixXd>(pca_res["h"]);
  MatrixXd eigf = as<MatrixXd>(pca_res["eig.functions"]);
  
  // Calculate S matrix
  MatrixXd S = H * (H.transpose() * H).inverse() * H.transpose();
  
  // Norm calculation
  VectorXd xnorm(n);
  for (int i = 0; i < n; i++) {
    xnorm(i) = sqrt(inpro(X.row(i), X.row(i), tm));
  }
  // Calculate mean of xnorm
  double mean_xnorm = xnorm.mean();
  
  // Calculate variance of xnorm
  double variance_xnorm = (xnorm.array() - mean_xnorm).square().sum() / (n - 1);
  
  // Calculate standard deviation
  double stddev_xnorm = sqrt(variance_xnorm);
  
  // Calculate h1
  double h1 = 2.34 * stddev_xnorm * pow(n, -0.2);

  // Distance matrix calculation
  MatrixXd d = MatrixXd::Constant(n, n, 0.0);
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      d(i, j) = sqrt(inpro(X.row(i) - X.row(j), X.row(i) - X.row(j), tm));
    }
  }
  // Kernel matrix
  MatrixXd Kz = MatrixXd::Constant(n, n, 0.0);
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      Kz(i, j) = fkernel(d(i, j), h1);
    }
  }
  Kz = Kz + Kz.transpose();
  
  // Estimate betahat
  MatrixXd I = MatrixXd::Identity(n, n);
  MatrixXd Zt = Z.transpose();
  VectorXd betahat = (Zt * (I - S) * Z).inverse() * Zt * (I - S) * Y;
  // Estimate alphahat
  MatrixXd alphahat = (H.transpose() * H).inverse() * H.transpose() * (Y - Z * betahat);
  VectorXd alphahatf = eigf * alphahat;
  
  // Estimate ehat
  VectorXd ehat(n);
  for (int i = 0; i < n; i++) {
    ehat(i) = Y(i) - Z.row(i) * betahat - inpro(alphahatf, X.row(i), tm);
  }
  // Rcpp::Rcout << Kz(3,1) << std::endl;
  // Final statistic calculation
  double statvalue = ehat.transpose() * Kz * ehat;
  return List::create(
    Named("stat") = statvalue,
    Named("ehat") = ehat, 
    Named("Kz") = Kz, 
    Named("alphahatf") = alphahatf, 
    Named("betahat") = betahat, 
    Named("H") = H, 
    Named("S") = S, 
    Named("eigf") = eigf
  );
}

// [[Rcpp::export]]
VectorXd BE_boot(
    const Eigen::MatrixXd &X, const Eigen::MatrixXd &Z, 
    const VectorXd &alphahatf, const VectorXd &betahat, 
    const MatrixXd &H, const MatrixXd &S, const MatrixXd &eigf, 
    const VectorXd &ehatcun, const MatrixXd &Kz, 
    int B, const VectorXd &tm
) {
  int n = Z.rows();
  VectorXd bvalue(B);
  
  for (int index = 0; index < B; index++) {
    VectorXd yb = VectorXd::Zero(n);
    VectorXd temp = VectorXd::Random(n);  // Uniform random numbers between -1 and 1
    //Rcpp::Rcout << temp << std::endl;
    temp = (temp.array() + 1.0) / 2.0;     // Normalize to [0, 1]
    
    VectorXd eb = (1 + sqrt(5)) / 2 - sqrt(5) * (temp.array() < (5 + sqrt(5)) / 10).cast<double>();
    VectorXd ebnew = eb.array() * ehatcun.array();
    
    for (int i = 0; i < n; i++) {
      yb(i) = Z.row(i) * betahat + inpro(alphahatf, X.row(i), tm) + ebnew(i);
    }
    
    VectorXd bbetahat = (Z.transpose() * (MatrixXd::Identity(n, n) - S) * Z).inverse() * (Z.transpose() * (MatrixXd::Identity(n, n) - S) * yb);
    VectorXd balphahat = (H.transpose() * H).inverse() * (H.transpose() * (yb - Z * bbetahat));
    VectorXd balphahatf = eigf * balphahat;
    
    VectorXd ehat = VectorXd::Zero(n);
    for (int i = 0; i < n; i++) {
      ehat(i) = yb(i) - Z.row(i) * bbetahat - inpro(balphahatf, X.row(i), tm);
    }
    
    bvalue(index) = ehat.transpose() * Kz * ehat;
  }
  
  return bvalue;
}


// [[Rcpp::export]]
List OUR_stat(
    const Eigen::MatrixXd &X, 
    const Eigen::VectorXd &Y, 
    const Eigen::MatrixXd &Z, 
    int m, 
    const Eigen::VectorXd &tm, 
    const Eigen::MatrixXd &covhat
) {
  int n = X.rows();
  
  // PCA results
  List pca_res = pca_rcppeigen(X, m, tm);
  MatrixXd H = as<MatrixXd>(pca_res["h"]);
  MatrixXd eigf = as<MatrixXd>(pca_res["eig.functions"]);
  
  // Calculate S matrix
  MatrixXd S = H * (H.transpose() * H).inverse() * H.transpose();
  
  // Norm calculation
  VectorXd xnorm(n);
  for (int i = 0; i < n; i++) {
    xnorm(i) = sqrt(inpro(X.row(i), X.row(i), tm));
  }
  // Calculate mean of xnorm
  double mean_xnorm = xnorm.mean();
  
  // Calculate variance of xnorm
  double variance_xnorm = (xnorm.array() - mean_xnorm).square().sum() / (n - 1);
  
  // Calculate standard deviation
  double stddev_xnorm = sqrt(variance_xnorm);
  
  // Calculate h1
  double h1 = 2.34 * stddev_xnorm * pow(n, -0.2);
  
  // Distance matrix calculation
  MatrixXd d = MatrixXd::Constant(n, n, 0.0);
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      d(i, j) = sqrt(inpro(X.row(i) - X.row(j), X.row(i) - X.row(j), tm));
    }
  }
  // Kernel matrix
  MatrixXd Kz = MatrixXd::Constant(n, n, 0.0);
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      Kz(i, j) = fkernel(d(i, j), h1);
    }
  }
  Kz = Kz + Kz.transpose();
  
  // Estimate betahat
  MatrixXd I = MatrixXd::Identity(n, n);
  MatrixXd Zt = Z.transpose();
  VectorXd betahat = (Zt * (I - S) * Z - n * covhat).inverse() * Zt * (I - S) * Y;
  // Estimate alphahat
  MatrixXd alphahat = (H.transpose() * H).inverse() * H.transpose() * (Y - Z * betahat);
  VectorXd alphahatf = eigf * alphahat;
  
  // Estimate ehat
  VectorXd ehat(n);
  for (int i = 0; i < n; i++) {
    ehat(i) = Y(i) - Z.row(i) * betahat - inpro(alphahatf, X.row(i), tm);
  }
  // Rcpp::Rcout << Kz(3,1) << std::endl;
  // Final statistic calculation
  double statvalue = ehat.transpose() * Kz * ehat;
  return List::create(
    Named("stat") = statvalue,
    Named("ehat") = ehat, 
    Named("Kz") = Kz, 
    Named("alphahatf") = alphahatf, 
    Named("betahat") = betahat, 
    Named("H") = H, 
    Named("S") = S, 
    Named("eigf") = eigf
  );
}

// [[Rcpp::export]]
VectorXd OUR_boot(
    const Eigen::MatrixXd &X, const Eigen::MatrixXd &Z, 
    const VectorXd &alphahatf, const VectorXd &betahat, 
    const MatrixXd &H, const MatrixXd &S, const MatrixXd &eigf, 
    const VectorXd &ehatcun, const MatrixXd &Kz, 
    int B, const VectorXd &tm
) {
  int n = Z.rows();
  VectorXd bvalue(B);
  
  for (int index = 0; index < B; index++) {
    VectorXd yb = VectorXd::Zero(n);
    VectorXd temp = VectorXd::Random(n);  // Uniform random numbers between -1 and 1
    //Rcpp::Rcout << temp << std::endl;
    temp = (temp.array() + 1.0) / 2.0;     // Normalize to [0, 1]
    
    VectorXd eb = (1 + sqrt(5)) / 2 - sqrt(5) * (temp.array() < (5 + sqrt(5)) / 10).cast<double>();
    VectorXd ebnew = eb.array() * ehatcun.array();
    
    for (int i = 0; i < n; i++) {
      yb(i) = Z.row(i) * betahat + inpro(alphahatf, X.row(i), tm) + ebnew(i);
    }
    
    VectorXd bbetahat = (Z.transpose() * (MatrixXd::Identity(n, n) - S) * Z).inverse() * (Z.transpose() * (MatrixXd::Identity(n, n) - S) * yb);
    VectorXd balphahat = (H.transpose() * H).inverse() * (H.transpose() * (yb - Z * bbetahat));
    VectorXd balphahatf = eigf * balphahat;
    
    VectorXd ehat = VectorXd::Zero(n);
    for (int i = 0; i < n; i++) {
      ehat(i) = yb(i) - Z.row(i) * bbetahat - inpro(balphahatf, X.row(i), tm);
    }
    
    bvalue(index) = ehat.transpose() * Kz * ehat;
  }
  
  return bvalue;
}

// Function to remove the j-th row from a matrix
MatrixXd removeRow(const MatrixXd &matrix, int j) {
  MatrixXd result(matrix.rows() - 1, matrix.cols());
  result.topRows(j) = matrix.topRows(j);
  result.bottomRows(matrix.rows() - j - 1) = matrix.bottomRows(matrix.rows() - j - 1);
  return result;
}

// Function to remove the j-th element from a vector
VectorXd removeElement(const VectorXd &vector, int j) {
  VectorXd result(vector.size() - 1);
  result.head(j) = vector.head(j);
  result.tail(vector.size() - j - 1) = vector.tail(vector.size() - j - 1);
  return result;
}

// Main function to perform cross-validation
// [[Rcpp::export]]
MatrixXd cross_validation(
    int rep, int n, double cov_u, double delta, const VectorXd &tm
) {
  MatrixXd cvresult = MatrixXd::Zero(rep, 10);
  
  // Define R function `DGP`
  Rcpp::Function DGP("DGP");
  MatrixXd diaguu = MatrixXd::Identity(2, 2) * cov_u;
  for (int i = 0; i < rep; i++) {
    Rcpp::Rcout << "loop: " << i+1 << "\r";
    // Call `generate` function in R and convert output to Eigen types
    List data = DGP(n, cov_u, delta);
    MatrixXd X = as<MatrixXd>(data["X"]);
    VectorXd Y = as<VectorXd>(data["Y"]);
    MatrixXd W = as<MatrixXd>(data["W"]);
    
    for (int m = 1; m <= 10; m++) {
      for (int j = 0; j < n; j++) {

        VectorXd x = X.row(j);
        double y = Y(j);
        VectorXd w = W.row(j);
        
        MatrixXd datax = removeRow(X, j);   // Exclude j-th row from X
        VectorXd datay = removeElement(Y, j);   // Exclude j-th element from Y
        MatrixXd dataw = removeRow(W, j);   // Exclude j-th row from W
        
        // PCA calculation using Rcpp function `pca_rcppeigen`
        List pca_res = pca_rcppeigen(datax, m, tm);
        MatrixXd H = as<MatrixXd>(pca_res["h"]);
        MatrixXd eigf = as<MatrixXd>(pca_res["eig.functions"]);
        
        MatrixXd S = H * ((H.transpose() * H).inverse()) * H.transpose();
        
        // Calculate betahat
        MatrixXd WTW = dataw.transpose() * (MatrixXd::Identity(n-1, n-1) - S) * dataw - (n - 1) * diaguu;
        VectorXd WTy = dataw.transpose() * (MatrixXd::Identity(n-1, n-1) - S) * datay;
        VectorXd betahat = WTW.inverse() * WTy;
        
        // Calculate alphahat and alphahatf
        VectorXd alphahat = (H.transpose() * H).inverse() * H.transpose() * (datay - dataw * betahat);
        VectorXd alphahatf = eigf * alphahat;
        
        // Calculate prediction error and update cvresult
        double prediction_error = pow(y - w.dot(betahat) - inpro(alphahatf, x, tm), 2) - betahat.transpose() * diaguu * betahat;
        cvresult(i, m-1) += prediction_error;
      }
    }
  }
  
  return cvresult;
}





