#include "R_ext/RS.h" /* for e.g., `F77_CALL` macro */
#include "RcppArmadillo.h"
#include "BLAS_LINPACK.h"
#include <memory>

extern "C" {
  void F77_NAME(dchud)(
      double*, int*, int*, double*, double*, int*, int*, double*, double*,
      double*, double*);

  void F77_NAME(dchdd)(
      double*, int*, int*, double*, double*, int*, int*, double*, double*,
      double*, double*, int*);
}

//' @import Rcpp
//' @useDynLib rollRegres, .registration = TRUE
// [[Rcpp::export]]
arma::mat roll_cpp(const arma::mat &X, const arma::vec &Y, int window){
  int n = X.n_rows, p = X.n_cols;
  const int p_cnst = p;
  arma::mat X_T = X.t();

  // initalize output
  arma::mat out(p, n); // notice other order
  std::fill(out.begin(), out.end(), NA_REAL);

  // compute coefs
  bool is_first = true;
  double d_one = 1, ddum;
  int i_one = 1L, i_zero = 0L;
  std::unique_ptr<int    []> jpvt (new int[p]            );
  std::unique_ptr<double []> qraux(new double[p]         );
  std::unique_ptr<double []> X_qr (new double[window * p]);
  std::unique_ptr<double []> XtY  (new double[p]         );
  std::unique_ptr<double []> c    (new double[p]         );
  std::unique_ptr<double []> s    (new double[p]         );
  double *X_T_begin = X_T.begin();

  for(int i = 0; i < p; ++i){
    jpvt[i] = 0;
    XtY[i] = 0;
  }

  for(int i = window - 1; i < n; ++i){
    if(is_first){ // setup QR decomposition
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
      dqrdc(&X_qr[0], &window, &window, &p, &qraux[0], &jpvt[0], &work, &job);

    } else {
      // update R
      int inc_new = i * p, inc_old = (i - window) * p;

      F77_CALL(dchud)(
          &X_qr[0], &window, &p, X_T_begin + inc_new,
          &ddum, &i_zero, &i_zero, &ddum, &ddum, &c[0], &s[0]);

      int info;
      F77_CALL(dchdd)(
          &X_qr[0], &window, &p, X_T_begin + inc_old,
          &ddum, &i_zero, &i_zero, &ddum, &ddum, &c[0], &s[0], &info);

      if(info != 0)
        Rcpp::stop("'dchdd' failed");

      // update X^T y
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
    dtrsm(
        "L", "U", "T", "N", &p_cnst, &i_one /* Y has one column */, &d_one,
        &X_qr[0], &window /* LDA */,
        out.begin() + i * p, &p_cnst);
    dtrsm(
        "L", "U", "N" /*  only difference */, "N", &p_cnst, &i_one, &d_one,
        &X_qr[0], &window,
        out.begin() + i * p, &p_cnst);

    /* norm = 0;
    for(int k = i - (window - 1L); k <= i; ++k){
      double res = Y[k];
      for(int j = 0; j < p; ++j)
        res -= out[i * p + j] * X_T[k * p + j];
      norm += res * res;
    }
    norm = 1000; /* norm = std::sqrt(norm); */
  }

  return out.t();
}
