#include <R_ext/RS.h> /* for F77_CALL */
#include <R_ext/Linpack.h>
#include "RcppArmadillo.h"

//' @import Rcpp
//' @useDynLib rollRegres, .registration = TRUE
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix roll_cpp(
    const arma::vec &Y, const arma::mat &X,
    int window){
  int n = X.n_rows, p = X.n_cols;
  const int p_cnst = p;
  arma::mat X_T = X.t();

  // initalize output
  Rcpp::NumericMatrix out(p, n); // notice other order
  std::fill(out.begin(), out.end(), NA_REAL);

  // compute coefs
  bool is_first = true;
  const int dum_int = 1;
  const double dum_dub = 1;
  int    jpvt [p];
  double qraux[p];
  double X_qr [window * p];
  double XtY  [p];
  double c    [p];
  double dum  [window];

  for(int i = 0; i < p; ++i){
    jpvt[i] = 0;
    XtY[i] = 0;
  }

  for(int i = window - 1; i < n; ++i){
    if(is_first){ // setup QR decomposition
      // see example in 9.10 of Linpack user guide
      is_first = false;
      int job = 0; /* don't pivot */
      double work; /* not referenced when not pivoting */

      // copy values
      for(int j = 0; j < p; ++j){ /* for each covaraite */
        for(int k = 0; k < window; ++k){ /* for each observation */
          X_qr[k + j * window] = X[k + j * n];
          XtY[j] += Y[k] * X[k + j * n];
        }
      }

      /* compute qR and get get upper triangular matrix, R, in the QR
       * decomposition see `qr.default`, `qr.R`, `r-source/src/appl/dqrdc2.f`
       * We use `dqrdc` instead where we can select not to perform pivoting.
       * NOTICE: no rank check
       */
      F77_CALL(dqrdc)(
          &X_qr[0], &window, &window, &p, &qraux[0], &jpvt[0], &work, &job);

      Rcpp::Rcout << 1 << std::endl;

    } else {
      // update R
      int inc_new = i * p, inc_old = (i - window) * p;

      Rcpp::Rcout << 2 << " " << inc_new << " " << inc_old << std::endl;
      F77_CALL(dchud)(
          &X_qr[0], &window, &p, X_T.begin() + inc_new,
          &dum[0], 0L, 0L, &dum[0], &dum[0], &c[0], &dum[0]);

      Rcpp::Rcout << 3 << std::endl;
      int info;
      F77_CALL(dchdd)(
          &X_qr[0], &window, &p, X_T.begin() + inc_old,
          &dum[0], 0L, 0L, &dum[0], &dum[0], &c[0], &dum[0], &info);

      if(info != 0)
        Rcpp::stop("'dchdd' failed");

      // update X^T y
      Rcpp::Rcout << 4 << std::endl;
      for(int j = 0; j < p; ++j)
        XtY[j] += Y[i] * X_T[j + inc_new] - Y[i - window] * X_T[j + inc_old];
    }

    // compute X^-T X = X^T y
    for(int j = 0; j < p; ++j)
      out[j + i * p] = XtY[j];
    /* See `r-source/src/library/base/R/backsolve.R` and
     * `r-source/src/main/array.c`
     *
     * First compute
     *     T(R) X = B
     *
     * then compute
     *       R Z = X
     */
    Rcpp::Rcout << 5 << std::endl;
    F77_CALL(dtrsm)(
        "L", "U", "T", "N", &p_cnst, &dum_int /* Y has one column */, &dum_dub,
        &X_qr[0], &window /* LDA */,
        out.begin() + i * p, &p_cnst);
    F77_CALL(dtrsm)(
        "L", "U", "N" /*  only difference */, "N", &p_cnst, &dum_int, &dum_dub,
        &X_qr[0], &window,
        out.begin() + i * p, &p_cnst);
  }

  return out;
}
