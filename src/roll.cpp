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

/* see https://stats.stackexchange.com/a/72215/81865 */
inline double
  update_sse(const double x, double &sse, int &t, double &x_bar){
    t += 1L;
    double e_t = x - x_bar;
    x_bar += e_t / t;
    sse += e_t * (x - x_bar);
    return sse;
}

/* TODO: assume that the reverse is also stable */
inline double
  downdate_sse(const double x, double &sse, int &t, double &x_bar){
    double e_t = x - x_bar;
    x_bar = (t / (t - 1.)) * x_bar -  x / (t - 1.);
    sse -= e_t * (x - x_bar);
    t -= 1L;
    return sse;
  }

inline double dot(const double *x, const double *y, const unsigned int p){
  double out = 0;
  for(unsigned int j = 0; j < p; ++j, ++x, ++y)
    out += *x * *y;
  return out;
}

//' @import Rcpp
//' @useDynLib rollRegres, .registration = TRUE
// [[Rcpp::export]]
Rcpp::List roll_cpp(
    const arma::mat &X, const arma::vec &Y, int window,
      const bool do_compute_R_sqs, const bool do_compute_sigmas,
      const bool do_1_step_forecasts){
  int n = X.n_rows, p = X.n_cols;
  const int p_cnst = p;
  arma::mat X_T = X.t();

  /* initalize output */
  arma::mat out(p, n); // notice other order
  arma::vec R_sqs, sigmas, one_step_forecasts;
  std::fill(out.begin()                 , out.end()               , NA_REAL);
  if(do_compute_R_sqs){
    R_sqs.set_size(n);
    std::fill(R_sqs.begin()             , R_sqs.end()             , NA_REAL);
  }
  if(do_compute_sigmas){
    sigmas.set_size(n);
    std::fill(sigmas.begin()            , sigmas.end()            , NA_REAL);
  }
  if(do_1_step_forecasts){
    one_step_forecasts.set_size(n);
    std::fill(one_step_forecasts.begin(), one_step_forecasts.end(), NA_REAL);
  }

  /* define intermediates */
  bool is_first = true;
  double d_one = 1, ddum, y_bar = 0., ss_tot = 0.;
  int i_one = 1L, i_zero = 0L, t = 0L;
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

  /* compute values */
  for(int i = window - 1; i < n; ++i){
    if(is_first){ // setup QR decomposition
      is_first = false;
      int job = 0; /* do not pivot */
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

      if(do_compute_R_sqs)
        while(t < window)
          update_sse(Y[t], ss_tot, t, y_bar);

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

      // update and downdate X^T y
      for(int j = 0; j < p; ++j)
        XtY[j] += Y[i] * X_T[j + inc_new] - Y[i - window] * X_T[j + inc_old];

      if(do_compute_R_sqs){
        update_sse  (Y[i         ], ss_tot, t, y_bar);
        downdate_sse(Y[i - window], ss_tot, t, y_bar);
      }
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

    /* compute other values if requested */
    double ss_reg = 0;
    if(do_compute_R_sqs or do_compute_sigmas){
      for(int k = i - (window - 1L); k <= i; ++k){
        double res = Y[k] - dot(&out[i * p], &X_T[k * p], p);
        ss_reg += res * res;
      }

      if(do_compute_sigmas)
        sigmas[i] = std::sqrt(ss_reg / (window - p));
    }

    if(do_compute_R_sqs){
      R_sqs[i] = (ss_tot  - ss_reg) / ss_tot;
    }

    if(do_1_step_forecasts and i < n - 1L){
      int next_i = i + 1L;
      one_step_forecasts[next_i] = dot(&out[i * p], &X_T[next_i * p], p);
    }
  }

  Rcpp::List out_list;
  out_list["coefs"] = out.t();
  if(do_compute_sigmas)
    out_list["sigmas"] = sigmas;
  else
    out_list["sigmas"] = R_NilValue;

  if(do_compute_R_sqs)
    out_list["r.squareds"] = R_sqs;
  else
    out_list["r.squareds"] = R_NilValue;

  if(do_1_step_forecasts)
    out_list["one_step_forecasts"] = one_step_forecasts;
  else
    out_list["one_step_forecasts"] = R_NilValue;

  return out_list;
}
