#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix beta_tau_mt(int n2, int m, int d,float delta_tau, NumericMatrix beta0_tau_hat,NumericVector tau_seq,
                          NumericMatrix u_mt, NumericMatrix int_u_mt, NumericMatrix Xmis_mt, NumericMatrix ymis_mt)
{
  int k = 0;
  float beta = 0;
  for(int l = 0; l < d; l++){
    for (int i = 0; i < n2; i++) {
      for (int j = 0; j < m; j++){
        k = int_u_mt(i,j);
        // 插值法确定beta
        beta = ( beta0_tau_hat(l ,k-1)*(tau_seq[k]- u_mt(i,j)) + beta0_tau_hat(l,k)*
          (u_mt(i,j) - tau_seq[k-1]))/delta_tau;
        // y_tilde
        ymis_mt( i, j) = ymis_mt( i, j) + Xmis_mt(i , l)*beta;
      }
    }
  }
  return(ymis_mt);
}
