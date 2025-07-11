// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// outerminus
arma::mat outerminus(arma::vec X);
RcppExport SEXP _rRAT_outerminus(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(outerminus(X));
    return rcpp_result_gen;
END_RCPP
}
// bubbleSortX
arma::mat bubbleSortX(arma::vec x, arma::vec y);
RcppExport SEXP _rRAT_bubbleSortX(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(bubbleSortX(x, y));
    return rcpp_result_gen;
END_RCPP
}
// mtau
double mtau(arma::vec x, arma::vec y);
RcppExport SEXP _rRAT_mtau(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(mtau(x, y));
    return rcpp_result_gen;
END_RCPP
}
// taucpp
arma::vec taucpp(arma::mat X);
RcppExport SEXP _rRAT_taucpp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(taucpp(X));
    return rcpp_result_gen;
END_RCPP
}
// allmain
arma::vec allmain(arma::mat X);
RcppExport SEXP _rRAT_allmain(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(allmain(X));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rRAT_outerminus", (DL_FUNC) &_rRAT_outerminus, 1},
    {"_rRAT_bubbleSortX", (DL_FUNC) &_rRAT_bubbleSortX, 2},
    {"_rRAT_mtau", (DL_FUNC) &_rRAT_mtau, 2},
    {"_rRAT_taucpp", (DL_FUNC) &_rRAT_taucpp, 1},
    {"_rRAT_allmain", (DL_FUNC) &_rRAT_allmain, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_rRAT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
