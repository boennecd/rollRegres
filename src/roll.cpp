#include "R_ext/RS.h" /* for e.g., `F77_CALL` macro */
#include "RcppArmadillo.h"
#include "BLAS_LINPACK.h"
#include <memory>
#include <cstring>
#include <sstream>

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
    /* TODO: large number minus small number -- likely not good */
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

inline int find_stard_end(
    int &start, int &end, const int delta_grp, const int *grp, const int n,
    /* additional arguments if stopping can also be caused by a sufficient
     * number of observations */
    const bool use_min_obs = false, const int min_obs = 0L, int nobs = 1L){
  if(start == n)
    Rcpp::stop("Invalid 'start' and 'n'");

  end = start;
  if(delta_grp <= 0L){
    return 0L;
  }

  int start_grp = *(grp + start), next_grp = start_grp, old_grp, length = 0L;

  bool do_stop = false;
  do {
    ++length; ++nobs;

    if(end == n - 1L){ /* reached the end */
      ++end;
      break;
    }

    old_grp = next_grp;
    next_grp = *(grp + ++end);
    if(old_grp > next_grp)
      Rcpp::stop("'grp' is not sorted");

    do_stop = next_grp - start_grp >= delta_grp;
    if(use_min_obs && next_grp != old_grp && nobs >= min_obs)
      do_stop = true;
  } while(!do_stop);

  return length;
}

//' @import Rcpp
//' @useDynLib rollRegres, .registration = TRUE
// [[Rcpp::export]]
Rcpp::List roll_cpp(
    const arma::mat &X, const arma::vec &Y, const int window,
    const bool do_compute_R_sqs,
    const bool do_compute_sigmas, const bool do_1_step_forecasts,
    arma::ivec grp, const bool use_grp, const bool do_downdates,
    const bool use_min_obs = false, const int min_obs = 0L){
  int n = X.n_rows, p = X.n_cols;
  const int p_cnst = p;
  arma::mat X_T = X.t();
  const double too_low = p * 3;

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
  bool is_first = true, do_warn = false;
  double d_one = 1, ddum, y_bar = 0., ss_tot = 0.;
  int i_one = 1L, i_zero = 0L, t = 0L;
  std::unique_ptr<int    []> jpvt (new int[p]   );
  std::unique_ptr<double []> qraux(new double[p]);
  std::unique_ptr<double []> X_qr;
  int ld_X_qr;
  std::unique_ptr<double []> XtY  (new double[p]);
  std::unique_ptr<double []> c    (new double[p]);
  std::unique_ptr<double []> s    (new double[p]);
  double *X_T_begin = X_T.begin();

  for(int i = 0; i < p; ++i){
    jpvt[i] = 0;
    XtY[i] = 0;
  }
  int start = 0, end = 0L,
    delete_start = (do_downdates) ? -1L : 0L, delete_end = 0L,
    sample_size = 0L, this_grp_start = -1L;

  /* compute values */
  for(int i = 0L; i < n; i = end){
    if(is_first){ // setup QR decomposition
      if(use_grp){
        sample_size += find_stard_end(
          start         , end, window, grp.begin(), n, use_min_obs, min_obs);
        /* need to find `this_grp_start` */
        int last_grp = grp[end - 1L];
        this_grp_start = end - 1L;
        while(grp[this_grp_start - 1L] == last_grp and this_grp_start > 0L)
          --this_grp_start;

      } else {
        start = 0L;
        end = window;
        this_grp_start = end - 1L;
        sample_size = window;

      }

      is_first = false;
      int job = 0; /* do not pivot */
      double work; /* not referenced when not pivoting */

      // copy values
      ld_X_qr = end - start;
      X_qr.reset(new double[p * ld_X_qr]);
      for(int j = 0; j < p; ++j){ /* for each covaraite */
        for(int k = start; k < end; ++k){ /* for each observation */
          X_qr[k + j * ld_X_qr] = X[k + j * n];
          XtY[j] += Y[k] * X[k + j * n];
        }
      }

      /* compute qR and get get upper triangular matrix, R, in the QR
       * decomposition see `qr.default`, `qr.R`, `r-source/src/appl/dqrdc2.f`
       * We use `dqrdc` instead where we can select not to perform pivoting.
       * NOTICE: no rank check
       */
      dqrdc(&X_qr[0], &ld_X_qr, &ld_X_qr, &p, &qraux[0], &jpvt[0],
            &work, &job);

      if(do_compute_R_sqs)
        while(t < end - start)
          update_sse(Y[t], ss_tot, t, y_bar);

    } else {
      /* find data start and end for both new and old data */
      if(use_grp){
        start = end;
        this_grp_start = start;
        sample_size += find_stard_end(start, end, 1L , grp.begin(), n);
        /* E.g., window 10, grp[start] is 15, grp[delete_start] is 3 so we
         * grp[delete_end] to be at least 6 as we want at most 6-15 so diff
         * should be between grp[delete_start] and grp[delete_end] should be
         *    3 = 15 - 3 - 10
         * or same example with 10, 11, and 0 which yields 1 = 11 - 0 - 10
         */
        if(do_downdates){
          delete_start = delete_end;
          sample_size -= find_stard_end(
            delete_start, delete_end,
            grp[start] - grp[delete_start] - window + 1L,
            grp.begin(), n);
        }

      } else {
        start = end;
        ++end;
        ++this_grp_start;
        if(do_downdates){
          ++delete_start;
          ++delete_end;

        } else
          ++sample_size;

      }

      /* updates */
      for(int k = start; k < end; ++k){
        int inc_new = k * p;
        F77_CALL(dchud)(
            &X_qr[0], &ld_X_qr, &p, X_T_begin + inc_new,
            &ddum, &i_zero, &i_zero, &ddum, &ddum, &c[0], &s[0]);

        /* update X^T y */
        for(int j = 0; j < p; ++j)
          XtY[j] += Y[k] * X_T[j + inc_new];

        if(do_compute_R_sqs)
          update_sse(Y[k], ss_tot, t, y_bar);
      }

      /* downdates */
      for(int k = delete_start; k < delete_end; ++k){
        int inc_old = k * p;
        int info;
        F77_CALL(dchdd)(
            &X_qr[0], &ld_X_qr, &p, X_T_begin + inc_old,
            &ddum, &i_zero, &i_zero, &ddum, &ddum, &c[0], &s[0], &info);

        if(info != 0){
          std::ostringstream os;
          os << "'dchdd' failed with code " << info;
          Rcpp::stop(os.str());
        }

        // downdate X^T y
        for(int j = 0; j < p; ++j)
          XtY[j] -= Y[k] * X_T[j + inc_old];

        if(do_compute_R_sqs)
          downdate_sse(Y[k], ss_tot, t, y_bar);
      }
    }

    if(!do_warn and sample_size < too_low)
      do_warn = true;

    // compute X^-T X = X^T y
    int coef_start = this_grp_start * p;
    for(int j = 0; j < p; ++j)
      out[j + coef_start] = XtY[j];

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
        &X_qr[0], &ld_X_qr /* LDA */,
        out.begin() + coef_start, &p_cnst);
    dtrsm(
        "L", "U", "N" /*  only difference */, "N", &p_cnst, &i_one, &d_one,
        &X_qr[0], &ld_X_qr,
        out.begin() + coef_start, &p_cnst);

    /* copy values */
    {
      double *o = out.begin(), *o_copy = o + coef_start;
      for(int k = this_grp_start; k < end; ++k)
        std::memcpy(o + k * p, o_copy, p * sizeof(double));
    }

    /* compute other values if requested */
    double ss_reg = 0;
    if(do_compute_R_sqs or do_compute_sigmas){
      for(int k = delete_end; k < end; ++k){
        double res = Y[k] - dot(&out[coef_start], &X_T[k * p], p);
        ss_reg += res * res;
      }

      if(do_compute_sigmas){
        double new_sigma = std::sqrt(ss_reg / (sample_size - p));
        for(int k = this_grp_start; k < end; ++k)
          sigmas[k] = new_sigma;
      }
    }

    if(do_compute_R_sqs){
      double new_R_sqs = (ss_tot  - ss_reg) / ss_tot;
      for(int k = this_grp_start; k < end; ++k)
        R_sqs[k] = new_R_sqs;
    }

    if(do_1_step_forecasts and end < n){
      int next_i = end;
      int next_grp = grp[next_i];
      for(; next_i < n and grp[next_i] == next_grp; ++next_i)
        one_step_forecasts[next_i] =
          dot(&out[coef_start], &X_T[next_i * p], p);
    }
  }

  if(do_warn){
    /* see https://stackoverflow.com/a/24557819/5861244 and
     * github.com/boennecd/rollRegres/issues/2 */
    Rcpp::Function warning("warning");
    warning("low sample size relative to number of parameters");
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


// [[Rcpp::export(name = ".find_chunks")]]
Rcpp::List chunk(const arma::ivec grp, const unsigned int width,
                 const unsigned int min_obs){
  /* Idea: we make one pass through `grp`. We keep track of the number of
   *       observations, current chunk number, the current starting index, and
   *       whether we are about to set a new chunk */
  bool is_new_chunk = true;
  unsigned int nobs = 1L, chunk_n = 0, istart = 0, grp_start = 0;
  const arma::uword n = grp.n_elem;
  std::vector<int> grp_idx_start, grp_idx_stop, has_value_start;
  grp_idx_start  .reserve(n);
  grp_idx_stop   .reserve(n);
  has_value_start.reserve(n);

  const int *g = grp.begin();
  /* hopefully not a used group... */
  int cur_grp = std::numeric_limits<int>::min(), first_grp = *grp.begin();

  bool has_window_length = false;
  for(unsigned int i = 0; i < grp.n_elem; ++g, ++i, ++nobs){
    bool is_new_grp = *g != cur_grp;
    has_window_length = has_window_length or *g - first_grp >= (int)width - 1L;
    cur_grp = *g;
    if(is_new_grp){
      grp_start = i;

      /* need to update istart */
      const int *g1 = grp.begin() + istart;
      for(unsigned int j = istart; j <= i; ++j, ++g1){
        if(*g - *g1 < (int)width){
          istart = j;
          break;
        } else
          --nobs; /* one less observation in the window */
      }
    }

    bool is_last = i == grp.n_elem - 1L,
      is_last_in_group = is_last || *g != *(g + 1L);
    if(!is_last_in_group or !has_window_length)
      continue;

    if(is_new_chunk && nobs >= min_obs){
      /* we have a new chunk */
      is_new_chunk = false;

      /* we need to set the values from istart and where those with values
       * start */
      ++chunk_n;
      has_value_start.push_back(grp_start + 1L);
      grp_idx_start.push_back(istart + 1L);

    } else if(!is_new_chunk && nobs < min_obs){
      /* need to find a new chunk */
      is_new_chunk = true;
      grp_idx_stop.push_back(grp_start);

    }

    if(!is_new_chunk && is_last){
      grp_idx_stop.push_back(i + 1L);

    }

  }

  return Rcpp::List::create(
    Rcpp::Named("grp_idx_start")   = Rcpp::wrap(grp_idx_start),
    Rcpp::Named("grp_idx_stop")    = Rcpp::wrap(grp_idx_stop),
    Rcpp::Named("has_value_start") = Rcpp::wrap(has_value_start));
}

