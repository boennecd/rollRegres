// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dchud_wrap
void dchud_wrap(arma::mat& r, int ldr, int p, arma::vec& x, arma::mat& z, int ldz, int nz, double y, double& rho, arma::vec& c, arma::vec& s);
RcppExport SEXP _rollRegres_dchud_wrap(SEXP rSEXP, SEXP ldrSEXP, SEXP pSEXP, SEXP xSEXP, SEXP zSEXP, SEXP ldzSEXP, SEXP nzSEXP, SEXP ySEXP, SEXP rhoSEXP, SEXP cSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type ldr(ldrSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type ldz(ldzSEXP);
    Rcpp::traits::input_parameter< int >::type nz(nzSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type s(sSEXP);
    dchud_wrap(r, ldr, p, x, z, ldz, nz, y, rho, c, s);
    return R_NilValue;
END_RCPP
}
// dchdd_wrap
void dchdd_wrap(arma::mat& r, int ldr, int p, arma::vec& x, arma::mat& z, int ldz, int nz, double y, double& rho, arma::vec& c, arma::vec& s, int& info);
RcppExport SEXP _rollRegres_dchdd_wrap(SEXP rSEXP, SEXP ldrSEXP, SEXP pSEXP, SEXP xSEXP, SEXP zSEXP, SEXP ldzSEXP, SEXP nzSEXP, SEXP ySEXP, SEXP rhoSEXP, SEXP cSEXP, SEXP sSEXP, SEXP infoSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type ldr(ldrSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type ldz(ldzSEXP);
    Rcpp::traits::input_parameter< int >::type nz(nzSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type s(sSEXP);
    Rcpp::traits::input_parameter< int& >::type info(infoSEXP);
    dchdd_wrap(r, ldr, p, x, z, ldz, nz, y, rho, c, s, info);
    return R_NilValue;
END_RCPP
}
// roll_cpp
Rcpp::List roll_cpp(const arma::mat& X, const arma::vec& Y, const int window, const bool do_compute_R_sqs, const bool do_compute_sigmas, const bool do_1_step_forecasts, arma::ivec grp, const bool use_grp, const bool do_downdates, const bool use_min_obs, const int min_obs);
RcppExport SEXP _rollRegres_roll_cpp(SEXP XSEXP, SEXP YSEXP, SEXP windowSEXP, SEXP do_compute_R_sqsSEXP, SEXP do_compute_sigmasSEXP, SEXP do_1_step_forecastsSEXP, SEXP grpSEXP, SEXP use_grpSEXP, SEXP do_downdatesSEXP, SEXP use_min_obsSEXP, SEXP min_obsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type window(windowSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_compute_R_sqs(do_compute_R_sqsSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_compute_sigmas(do_compute_sigmasSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_1_step_forecasts(do_1_step_forecastsSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type grp(grpSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_grp(use_grpSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_downdates(do_downdatesSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_min_obs(use_min_obsSEXP);
    Rcpp::traits::input_parameter< const int >::type min_obs(min_obsSEXP);
    rcpp_result_gen = Rcpp::wrap(roll_cpp(X, Y, window, do_compute_R_sqs, do_compute_sigmas, do_1_step_forecasts, grp, use_grp, do_downdates, use_min_obs, min_obs));
    return rcpp_result_gen;
END_RCPP
}
// chunk
Rcpp::List chunk(const arma::ivec grp, const unsigned int width, const unsigned int min_obs);
RcppExport SEXP _rollRegres_chunk(SEXP grpSEXP, SEXP widthSEXP, SEXP min_obsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::ivec >::type grp(grpSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type width(widthSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type min_obs(min_obsSEXP);
    rcpp_result_gen = Rcpp::wrap(chunk(grp, width, min_obs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rollRegres_dchud_wrap", (DL_FUNC) &_rollRegres_dchud_wrap, 11},
    {"_rollRegres_dchdd_wrap", (DL_FUNC) &_rollRegres_dchdd_wrap, 12},
    {"_rollRegres_roll_cpp", (DL_FUNC) &_rollRegres_roll_cpp, 11},
    {"_rollRegres_chunk", (DL_FUNC) &_rollRegres_chunk, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_rollRegres(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
