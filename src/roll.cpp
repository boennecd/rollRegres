#include "R_ext/RS.h" /* for e.g., `F77_CALL` macro */
#include "arma.h"
#include "BLAS_LINPACK.h"
#include <memory>
#include <cstring>
#include <sstream>
#include <array>

static constexpr int i_one = 1L, i_zero = 0L;
static constexpr double d_one = 1;
static constexpr char C_L = 'L', C_U = 'U', C_T = 'T', C_N = 'N';

extern "C" {
  void F77_NAME(dchud)(
      double*, const int*, const int*, const double*, double*, const int*,
      const int*, const double*, double*,
      double*, double*);
  void F77_NAME(dchdd)(
      double*, const int*, const int*, const double*, double*, const int*,
      const int*, const double*, double*, double*, double*, int*);
  void F77_NAME(dqrsl)(
      double*, int*, int*, int*, double*, const double*, double*, double*,
      double*, double*, double*, const int*, int*);
  void F77_NAME(dqrdc)(
      double*, int*, int*, int*, double*, int*, double*, const int*);
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

class roll_cpp_worker_base {
protected:
  const double *last_coef = nullptr;

public:
  /* updates the present solution to include new observations */
  virtual void update
  (const std::size_t, const std::size_t) = 0;

  /* downdate the present solution to exclude some observations */
  virtual void downdate
    (const std::size_t, const std::size_t) = 0;

  /* sets the coefficients to the passed pointer. Assumes that sufficient memory
   * is allocated. A pointer to the memory will be stored and may be used in
   * other member functions. */
  virtual void set_coef(double*) = 0;

  /* returns the sum of squared errors. May use a stored pointer to the last
   * computed coefficients. The class may not need the passed indices */
  virtual double get_ss(const std::size_t, const std::size_t) = 0;

  /* predicts the mean of a new outcome. May use a stored pointer to the last
   * computed coefficients */
  virtual double predict(const std::size_t) = 0;

  virtual ~roll_cpp_worker_base() = default;
};

/* Does essentially as described in the LINPACK guide */
class roll_cpp_worker_linpack final : public roll_cpp_worker_base {
  using iptr = std::unique_ptr<int[]>;
  using dptr = std::unique_ptr<double[]>;

  const arma::mat &X, X_T = X.t();
  const arma::vec &Y;
  const int n = X.n_rows;
  /* only non-const due to linpack definition */
  int p = X.n_cols;

  dptr
    /* upper upper triangular matrix with cholesky decomposition */
    R     = dptr(new double[p * p]),
    /* contains the first p elements of Q^\top y */
    z     = dptr(new double[p]),
    /* needed in Fortran calls */
    qraux = dptr(new double[p]),
    s     = dptr(new double[p]),
    c     = dptr(new double[p]);

  /* square root of sum of squared errors */
  double ss_sqrt;

  void check_start_end
  (const std::size_t start, const std::size_t end)
  {
#ifdef ROLL_DEBUG
    if(start >= end or (int)end > n)
      throw std::invalid_argument("Invalid 'start' and/or 'stop'");
#endif
  }

public:
  roll_cpp_worker_linpack
  (const arma::mat &X, const arma::vec &Y , const std::size_t start,
   const std::size_t end):
  X(X), Y(Y)
  {
    check_start_end(start, end);

    /* compute the cholesky decomposition. First, we need to setup a few
     * objects */
    iptr jpvt  = iptr(new int[p]);
    for(int i = 0; i < p; ++i)
      jpvt[i] = 0;

    static constexpr int job = 0; /* do not pivot */
    double work; /* not referenced when not pivoting */

    // copy values
    int ld_X_qr = end - start;
    dptr X_qr(new double[p * ld_X_qr]);
    for(int j = 0; j < p; ++j)
      for(unsigned int k = start; k < end; ++k)
        X_qr[k + j * ld_X_qr] = X[k + j * n];

    /* compute QR and get get upper triangular matrix, R, in the QR
     * decomposition. See `qr.default`, `qr.R`, `r-source/src/appl/dqrdc2.f`
     * We use `dqrdc` instead where we can select not to perform pivoting.
     * NOTICE: no rank check
     */
    F77_CALL(dqrdc)(
        &X_qr[0], &ld_X_qr, &ld_X_qr, &p, &qraux[0], &jpvt[0],
        &work, &job);

    /* copy R */
    {
      double *r = R.get();
      for(int i = 0; i < p; ++i){
        const double *xqr = X_qr.get() + i * ld_X_qr;
        for(int j = 0; j < p; ++r, ++xqr, ++j)
          *r = *xqr;
      }
    }

    /* get first part of Q^\top y and the squared norm */
    dptr Q_t_Y(new double[ld_X_qr]);
    static constexpr int job_qrsl = 1000L;
    double dummy = 0.;
    int info;

    F77_CALL(dqrsl)(
      X_qr.get(), &ld_X_qr, &ld_X_qr, &p, qraux.get(), Y.memptr() + start,
      &dummy, Q_t_Y.get(), &dummy, &dummy, &dummy, &job_qrsl, &info);

    if(info != 0L)
      throw std::runtime_error(
          "'dqrsl' failed with code " + std::to_string(info));

    /* compute ss */
    ss_sqrt = 0.;
    for(int i = p; i < ld_X_qr; ++i)
      ss_sqrt += *(Q_t_Y.get() + i) * *(Q_t_Y.get() + i);
    ss_sqrt = std::sqrt(ss_sqrt);

    /* store first p elements of Q^\top y */
    {
      for(int i = 0; i < p; ++i)
        *(z.get() + i) = *(Q_t_Y.get() + i);
    }
  }

  void update
  (const std::size_t start, const std::size_t end) override
  {
    check_start_end(start, end);

    const double *X_t_p, *y_p;
    unsigned int k;

    for(k = start, X_t_p = X_T.begin() + start * p, y_p = Y.memptr() + start;
        k < end; ++k, X_t_p += p, ++y_p)
      F77_CALL(dchud)(
          R.get(), &p, &p, X_t_p,
          z.get(), &p, &i_one, y_p, &ss_sqrt, &c[0], &s[0]);
  }

  void downdate
  (const std::size_t start, const std::size_t end) override
  {
    check_start_end(start, end);

    const double *X_t_p, *y_p;
    unsigned int k;
    int info;

    for(k = start, X_t_p = X_T.begin() + start * p, y_p = Y.memptr() + start;
        k < end; ++k, X_t_p += p, ++y_p){
      F77_CALL(dchdd)(
          R.get(), &p, &p, X_t_p,
          z.get(), &p, &i_one, y_p, &ss_sqrt, &c[0], &s[0], &info);

      if(info != 0)
        throw std::runtime_error(
            "'dchdd' failed with code " + std::to_string(info));
    }
  }

  void set_coef(double *out) override {
    last_coef = out;

    /* copy first elements of Q^\top y */
    for(int j = 0; j < p; ++j)
      *(out + j) = z[j];

    dtrsm(
      &C_L, &C_U, &C_N, &C_N, &p, &i_one, &d_one, R.get(), &p,
      out, &p);
  }

  double get_ss(const std::size_t start, const std::size_t end) override {
    return ss_sqrt * ss_sqrt;
  }

  double predict(const std::size_t i_new) override {
#ifdef ROLL_DEBUG
    if(!last_coef)
      throw std::logic_error("'predict' called before 'set_coef'");
    if((int)i_new >= n)
      throw std::invalid_argument("Invalid 'i_new'");
#endif
    return dot(last_coef, X_T.memptr() + i_new * p, p);
  }
};

/* The first implementation I made which stores X^\top y instead of the first
 * p rows of Q^\top y where Q is from the QR decompsotion of the design matrix */
class roll_cpp_worker final : public roll_cpp_worker_base {
  using iptr = std::unique_ptr<int[]>;
  using dptr = std::unique_ptr<double[]>;

  const arma::mat &X, X_T = X.t();
  const arma::vec &Y;
  const int n = X.n_rows;
  /* only non-const due to linpack definition */
  int p = X.n_cols, ld_X_qr;

  iptr jpvt  = iptr(new int[p]);
  dptr qraux = dptr(new double[p]),
       XtY   = dptr(new double[p]),
       s     = dptr(new double[p]),
       c     = dptr(new double[p]),
       X_qr;

  void check_start_end
    (const std::size_t start, const std::size_t end)
  {
#ifdef ROLL_DEBUG
    if(start >= end or (int)end > n)
      throw std::invalid_argument("Invalid 'start' and/or 'stop'");
#endif
  }

public:
  roll_cpp_worker
  (const arma::mat &X, const arma::vec &Y , const std::size_t start,
   const std::size_t end):
  X(X), Y(Y), ld_X_qr(end - start)
  {
    check_start_end(start, end);

    /* compute the cholesky decomposition. First, we need to setup a few
     * objects */
    for(int i = 0; i < p; ++i){
      jpvt[i] = 0;
      XtY[i] = 0;
    }

    int job = 0; /* do not pivot */
    double work; /* not referenced when not pivoting */

    // copy values
    X_qr.reset(new double[p * ld_X_qr]);
    for(int j = 0; j < p; ++j){ /* for each covaraite */
     for(unsigned int k = start; k < end; ++k){ /* for each observation */
        X_qr[k + j * ld_X_qr] = X[k + j * n];
        XtY[j] += Y[k] * X[k + j * n];
      }
    }

    /* compute QR and get get upper triangular matrix, R, in the QR
     * decomposition. See `qr.default`, `qr.R`, `r-source/src/appl/dqrdc2.f`
     * We use `dqrdc` instead where we can select not to perform pivoting.
     * NOTICE: no rank check
     */
    F77_CALL(dqrdc)(
        &X_qr[0], &ld_X_qr, &ld_X_qr, &p, &qraux[0], &jpvt[0],
        &work, &job);
  }

  void update
      (const std::size_t start, const std::size_t end) override
  {
    check_start_end(start, end);

    double ddum;
    const double *X_t_p;
    unsigned int k;

    for(k = start, X_t_p = X_T.begin() + start * p;
        k < end; ++k, X_t_p += p){
      F77_CALL(dchud)(
          &X_qr[0], &ld_X_qr, &p, X_t_p,
          &ddum, &i_zero, &i_zero, &ddum, &ddum, &c[0], &s[0]);

      /* update X^T y */
      for(int j = 0; j < p; ++j)
        XtY[j] += Y[k] * *(X_t_p + j);
    }
  }

  void downdate
  (const std::size_t start, const std::size_t end) override
  {
    check_start_end(start, end);

    double ddum;
    const double *X_t_p;
    unsigned int k;
    int info;

    for(k = start, X_t_p = X_T.begin() + start * p; k < end; ++k, X_t_p += p){
      F77_CALL(dchdd)(
          &X_qr[0], &ld_X_qr, &p, X_t_p,
          &ddum, &i_zero, &i_zero, &ddum, &ddum, &c[0], &s[0], &info);

      if(info != 0)
        throw std::runtime_error(
            "'dchdd' failed with code " + std::to_string(info));

      // downdate X^T y
      for(int j = 0; j < p; ++j)
        XtY[j] -= Y[k] * *(X_t_p + j);
    }
  }

  void set_coef(double *out) override {
    last_coef = out;

    /* copy XtY */
    for(int j = 0; j < p; ++j)
      *(out + j) = XtY[j];

    /* solve */
    auto solve_func = [&](char transpose){
      dtrsm(
        &C_L, &C_U, &transpose, &C_N, &p, &i_one, &d_one, &X_qr[0], &ld_X_qr,
        out, &p);
    };

    solve_func(C_T);
    solve_func(C_N);
  }

  double get_ss(const std::size_t start, const std::size_t end) override {
#ifdef ROLL_DEBUG
    if(!last_coef)
      throw std::logic_error("'get_ss' called before 'set_coef'");
#endif

    double ss_reg = 0.;
    const double *X_t_p = X_T.memptr() + start * p;
    for(unsigned int k = start; k < end; ++k, X_t_p += p){
      const double res = Y[k] - dot(last_coef, X_t_p, p);
      ss_reg += res * res;
    }

    return ss_reg;
  }

  double predict(const std::size_t i_new) override {
#ifdef ROLL_DEBUG
    if(!last_coef)
      throw std::logic_error("'predict' called before 'set_coef'");
    if((int)i_new >= n)
      throw std::invalid_argument("Invalid 'i_new'");
#endif
    return dot(last_coef, X_T.memptr() + i_new * p, p);
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
  const std::size_t n = X.n_rows, p = X.n_cols, too_low = p * 3L;

  /* initalize output */
  auto set_num_vec = [&](const bool not_empty){
    if(!not_empty)
      return Rcpp::NumericVector();
    Rcpp::NumericVector out(n);
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  };
  Rcpp::NumericMatrix out(p, n); // notice other order
  std::fill(out.begin(), out.end(), NA_REAL);

  Rcpp::NumericVector
    R_sqs              = set_num_vec(do_compute_R_sqs),
    sigmas             = set_num_vec(do_compute_sigmas),
    one_step_forecasts = set_num_vec(do_1_step_forecasts);

  /* define intermediates */
  bool is_first = true, do_warn = false;
  double y_bar = 0., ss_tot = 0.;
  unsigned int t = 0L;

  roll_cpp_indices idxs = ([&]{
    const arma::ivec *ptr = use_grp ? &grp : nullptr;
    return roll_cpp_indices(n, window, ptr, use_min_obs, min_obs);
  })();

  std::size_t
    delete_start = do_downdates ? -1L : 0L, delete_end = 0L,
    sample_size = 0L, this_grp_start = -1L;

  std::unique_ptr<roll_cpp_worker_base> worker;

  /* compute values */
  for(unsigned int i = 0L; i < n; i = idxs.end()){
    if(is_first){ // setup worker
      if(use_grp){
        /* need to find `this_grp_start` */
        int last_grp = grp[idxs.end() - 1L];
        this_grp_start = idxs.end() - 1L;
        while(grp[this_grp_start - 1L] == last_grp and this_grp_start > 0L)
          --this_grp_start;

      } else
        this_grp_start = idxs.end() - 1L;
      sample_size = idxs.end() - idxs.start();

      is_first = false;

      worker.reset(new roll_cpp_worker_linpack(
          X, Y, idxs.start(), idxs.end()));

      if(do_compute_R_sqs)
        while(t < sample_size)
          update_sse(Y[t], ss_tot, t, y_bar);

    } else {
      /* find data start and end for both new and old data */
      auto old = idxs.move(do_downdates);
      this_grp_start = old.end;
      sample_size = idxs.end() - idxs.start();

      /* updates */
      worker->update(old.end, idxs.end());
      for(unsigned int k = old.end; k < idxs.end() and do_compute_R_sqs;
          ++k)
        update_sse(Y[k], ss_tot, t, y_bar);

      /* downdates */
      if(do_downdates){
        delete_start = old.start.idx;
        delete_end = idxs.start();

        if(delete_start < delete_end)
          worker->downdate(delete_start, delete_end);

        for(unsigned int k = delete_start;
            k < delete_end and do_compute_R_sqs; ++k)
          downdate_sse(Y[k], ss_tot, t, y_bar);
      }
    }

    if(!do_warn and sample_size < too_low)
      do_warn = true;

    int coef_start = this_grp_start * p;
    worker->set_coef(out.begin() + coef_start);

    /* copy values */
    {
      double *o = out.begin(), *o_copy = o + coef_start;
      for(unsigned int k = this_grp_start; k < idxs.end(); ++k)
        std::memcpy(o + k * p, o_copy, p * sizeof(double));
    }

    /* compute other values if requested */
    if(do_compute_R_sqs or do_compute_sigmas){
      double ss_reg = worker->get_ss(idxs.start(), idxs.end());

      if(do_compute_sigmas){
        const double new_sigma = std::sqrt(ss_reg / (sample_size - p));
        for(unsigned int k = this_grp_start; k < idxs.end(); ++k)
          sigmas[k] = new_sigma;
      }

      if(do_compute_R_sqs){
        const double new_R_sqs = (ss_tot  - ss_reg) / ss_tot;
        for(unsigned int k = this_grp_start; k < idxs.end(); ++k)
          R_sqs[k] = new_R_sqs;
      }
    }

    if(do_1_step_forecasts and idxs.end() < n){
      unsigned int next_i = idxs.end();
      int next_grp = grp[next_i];
      for(; next_i < n and grp[next_i] == next_grp; ++next_i)
        one_step_forecasts[next_i] = worker->predict(next_i);
    }
  }

  if(do_warn){
    /* see https://stackoverflow.com/a/24557819/5861244 and
     * github.com/boennecd/rollRegres/issues/2 */
    Rcpp::Function warning("warning");
    warning("low sample size relative to number of parameters");
  }

  Rcpp::List out_list;
  out = Rcpp::transpose(out);
  out_list["coefs"] = std::move(out);
  if(do_compute_sigmas)
    out_list["sigmas"] = std::move(sigmas);
  else
    out_list["sigmas"] = R_NilValue;

  if(do_compute_R_sqs)
    out_list["r.squareds"] = std::move(R_sqs);
  else
    out_list["r.squareds"] = R_NilValue;

  if(do_1_step_forecasts)
    out_list["one_step_forecasts"] = std::move(one_step_forecasts);
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
