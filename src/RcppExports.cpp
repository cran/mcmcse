// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

// mbmc
mat mbmc(const mat& chain, double b);
RcppExport SEXP mcmcse_mbmc(SEXP chainSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat >::type chain(chainSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    __result = Rcpp::wrap(mbmc(chain, b));
    return __result;
END_RCPP
}
// inseq
List inseq(mat M, bool adjust);
RcppExport SEXP mcmcse_inseq(SEXP MSEXP, SEXP adjustSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type M(MSEXP);
    Rcpp::traits::input_parameter< bool >::type adjust(adjustSEXP);
    rcpp_result_gen = Rcpp::wrap(inseq(M, adjust));
    return rcpp_result_gen;
END_RCPP
}
// msvec
mat msvec(const mat& chain, double b, String method);
RcppExport SEXP mcmcse_msvec(SEXP chainSEXP, SEXP bSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat >::type chain(chainSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    __result = Rcpp::wrap(msvec(chain, b, method));
    return __result;
END_RCPP
}
