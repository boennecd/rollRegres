#include "R_ext/RS.h" /* for e.g., `F77_CALL` macro */
#include "RcppArmadillo.h"
#include "BLAS_LINPACK.h"
#include <memory>
#include <cstring>
#include <sstream>
#include <array>

static constexpr int i_one = 1L, i_zero = 0L;
static constexpr double d_one = 1;

extern "C" {
  void F77_NAME(dchud)(
      double*, const int*, const int*, double*, double*, const int*,
      const int*, double*, double*,
      double*, double*);

  void F77_NAME(dchdd)(
      double*, const int*, const int*, double*, double*, const int*,
      const int*, double*, double*, double*, double*, int*);
}

/* see https://stats.stackexchange.com/a/72215/81865 */
inline double
  update_sse(const double x, double &sse, unsigned int &t, double &x_bar){
    t += 1L;
    double e_t = x - x_bar;
    x_bar += e_t / t;
    sse += e_t * (x - x_bar);
    return sse;
}

/* TODO: assume that the reverse is also stable */
inline double
  downdate_sse(const double x, double &sse, unsigned int &t, double &x_bar){
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

/* class to keep track of indices when adding or removing observations */
class roll_cpp_indices {
  /* number of observations and window size */
  const std::size_t n, window;
  /* optional vector of group id for each observation */
  const arma::ivec *grp;
  /* is group indices supplied */
  const bool has_grp = grp, use_min_obs;
  /* optional minimum number of observations in first window */
  const std::size_t min_obs;

  /* move end one step. Either in terms of rows or the group id. The first
   * index is handled as a sepcial case */
  void move_end(){
    const bool is_first_call = end_.idx < 1L;

    if(!has_grp){
      if(is_first_call)
        end_.idx = std::min(start_.idx + window, n);
      else
        ++end_.idx;

      return;
    }

#ifdef ROLL_DEBUG
    if(has_grp and (!start_.grp or !end_.grp))
      throw std::runtime_error("grp not set");
#endif

    auto &i = end_.idx;
    auto &e_grp = end_.grp;
    const bool use_min_obs_now =
      /* only done on the first call */
      use_min_obs and is_first_call;

    do {
      ++i;
      ++e_grp;
#ifdef ROLL_DEBUG
      if(i < n and *e_grp - *start_.grp < 0L)
        throw std::runtime_error("groups are not sorted");
#endif
    } while(
      /* not reached the end */
      i < n and
      /* the difference between the groups is too small */
      (is_first_call and *e_grp - *start_.grp  < (int)window) and
      /* we do not use the minimum number of observations as a criteria or
       * the minimum number of observations is not reached */
      (!use_min_obs_now or i - start_.idx < min_obs));

    /* may have to move the end index to the start of next group */
    while(i < n and *e_grp == *(e_grp - 1L)){
      ++i;
      ++e_grp;
    }
  }

public:
  /* struct that contains the current index and optional pointer to group
   * index */
  struct index {
    std::size_t idx;
    const int *grp;

    operator std::size_t() const {
      return idx;
    }

    index(std::size_t idx, int *grp): idx(idx), grp(grp) { }
    index(std::size_t idx): index(idx, nullptr) { }
  };

private:
  /* current start and end index */
  index start_, end_;

public:
  const index& start() const {
    return start_;
  }
  const index& end() const {
    return end_;
  }

  roll_cpp_indices
    (const std::size_t n, const std::size_t window, const arma::ivec *grp,
     const bool use_min_obs, const std::size_t min_obs):
    n(n), window(window), grp(grp), use_min_obs(use_min_obs), min_obs(min_obs),
    start_(0L), end_(0L)
  {
#ifdef ROLL_DEBUG
    if(grp and grp->n_elem != n)
      throw std::invalid_argument("invalid 'grp' size and 'n'");
#endif

    if(grp)
      end_.grp = start_.grp = grp->begin();

    move_end();
  }

  struct old_idx_class {
    const index start, end;
  };

  /* move start and end one step and return copy of previous start and end */
  old_idx_class move(const bool do_move_start){
#ifdef ROLL_DEBUG
    if(end_.idx == n)
      throw std::logic_error("already at end");
#endif

    old_idx_class out { start_, end_ };

    if(do_move_start){
      if(has_grp){
        std::size_t &i   = start_.idx;
        auto &grp = start_.grp;
        const int next_grp_min = *end_.grp - window + 1L;

        while (
            /* have not reached the end*/
            i < n and
            /* group have not changed sufficiently */
            *grp < next_grp_min){
          ++i;
          ++grp;
        }

      } else
        ++start_.idx;

#ifdef ROLL_DEBUG
      if(start_.idx >= n)
        throw std::runtime_error("start is after end index");
#endif
    }

    move_end();

    return out;
  }
};

//' @import Rcpp
//' @useDynLib rollRegres, .registration = TRUE
// [[Rcpp::export]]
Rcpp::List roll_cpp(
    const arma::mat &X, const arma::vec &Y, const int window,
    const bool do_compute_R_sqs,
    const bool do_compute_sigmas, const bool do_1_step_forecasts,
    arma::ivec grp, const bool use_grp, const bool do_downdates,
    const bool use_min_obs = false, const int min_obs = 0L){
  const unsigned int n = X.n_rows;
  int p = X.n_cols;
  const int p_cnst = p;
  arma::mat X_T = X.t();
  const std::size_t too_low = p * 3L;

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
  double ddum, y_bar = 0., ss_tot = 0.;
  unsigned int t = 0L;
  std::unique_ptr<int    []> jpvt (new int[p]   );
  std::unique_ptr<double []> qraux(new double[p]);
  std::unique_ptr<double []> X_qr;
  int ld_X_qr;
  std::unique_ptr<double []> XtY  (new double[p]);
  std::unique_ptr<double []> c    (new double[p]);
  std::unique_ptr<double []> s    (new double[p]);
  double *X_T_begin = X_T.begin();

  roll_cpp_indices idxs = ([&]{
    const arma::ivec *ptr = use_grp ? &grp : nullptr;
    return roll_cpp_indices(n, window, ptr, use_min_obs, min_obs);
  })();

  for(int i = 0; i < p; ++i){
    jpvt[i] = 0;
    XtY[i] = 0;
  }
  std::size_t
    delete_start = do_downdates ? -1L : 0L, delete_end = 0L,
    sample_size = 0L, this_grp_start = -1L;

  /* compute values */
  for(unsigned int i = 0L; i < n; i = idxs.end()){
    if(is_first){ // setup QR decomposition
      if(use_grp){
        /* need to find `this_grp_start` */
        int last_grp = grp[idxs.end() - 1L];
        this_grp_start = idxs.end() - 1L;
        while(grp[this_grp_start - 1L] == last_grp and this_grp_start > 0L)
          --this_grp_start;

      } else
        this_grp_start = idxs.end() - 1L;

      is_first = false;
      int job = 0; /* do not pivot */
      double work; /* not referenced when not pivoting */

      // copy values
      sample_size = idxs.end() - idxs.start();
      ld_X_qr = sample_size;
      X_qr.reset(new double[p * ld_X_qr]);
      for(int j = 0; j < p; ++j){ /* for each covaraite */
        for(unsigned int k = idxs.start(); k < idxs.end(); ++k){ /* for each observation */
          X_qr[k + j * ld_X_qr] = X[k + j * n];
          XtY[j] += Y[k] * X[k + j * n];
        }
      }

      /* compute QR and get get upper triangular matrix, R, in the QR
       * decomposition see `qr.default`, `qr.R`, `r-source/src/appl/dqrdc2.f`
       * We use `dqrdc` instead where we can select not to perform pivoting.
       * NOTICE: no rank check
       */
      dqrdc(&X_qr[0], &ld_X_qr, &ld_X_qr, &p, &qraux[0], &jpvt[0],
            &work, &job);

      if(do_compute_R_sqs)
        while(t < sample_size)
          update_sse(Y[t], ss_tot, t, y_bar);

    } else {
      /* find data start and end for both new and old data */
      auto old = idxs.move(do_downdates);
      this_grp_start = old.end;
      sample_size = idxs.end() - idxs.start();

      if(do_downdates){
        delete_start = old.start.idx;
        delete_end = idxs.start();
      }

      /* updates */
      for(unsigned int k = old.end; k < idxs.end(); ++k){
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
      for(unsigned int k = delete_start; k < delete_end; ++k){
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
      for(unsigned int k = this_grp_start; k < idxs.end(); ++k)
        std::memcpy(o + k * p, o_copy, p * sizeof(double));
    }

    /* compute other values if requested */
    double ss_reg = 0;
    if(do_compute_R_sqs or do_compute_sigmas){
      for(unsigned int k = idxs.start(); k < idxs.end(); ++k){
        const double res = Y[k] - dot(&out[coef_start], &X_T[k * p], p);
        ss_reg += res * res;
      }

      if(do_compute_sigmas){
        const double new_sigma = std::sqrt(ss_reg / (sample_size - p));
        for(unsigned int k = this_grp_start; k < idxs.end(); ++k)
          sigmas[k] = new_sigma;
      }
    }

    if(do_compute_R_sqs){
      const double new_R_sqs = (ss_tot  - ss_reg) / ss_tot;
      for(unsigned int k = this_grp_start; k < idxs.end(); ++k)
        R_sqs[k] = new_R_sqs;
    }

    if(do_1_step_forecasts and idxs.end() < n){
      unsigned int next_i = idxs.end();
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

