/* export for testing and use in examples */
#include "R_ext/RS.h" /* for e.g., `F77_CALL` macro */
#include "RcppArmadillo.h"

extern "C" {
  void F77_NAME(dchud)(
      double*, int*, int*, double*, double*, int*, int*, double*, double*,
      double*, double*);

  void F77_NAME(dchdd)(
      double*, int*, int*, double*, double*, int*, int*, double*, double*,
      double*, double*, int*);
}

// [[Rcpp::export(rng = false)]]
void dchud_wrap
  (arma::mat &r, int ldr, int p, arma::vec &x, arma::mat &z, int ldz, int nz,
   double y, double &rho, arma::vec &c, arma::vec &s){
  F77_CALL(dchud)(
      r.begin(), &ldr, &p, x.begin(), z.begin(), &ldz, &nz, &y, &rho,
      c.begin(), s.begin());
}

// [[Rcpp::export(rng = false)]]
void dchdd_wrap
  (arma::mat &r, int ldr, int p, arma::vec &x, arma::mat &z, int ldz, int nz,
   double y, double &rho, arma::vec &c, arma::vec &s, int &info){
  F77_CALL(dchdd)(
      r.begin(), &ldr, &p, x.begin(), z.begin(), &ldz, &nz, &y, &rho,
      c.begin(), s.begin(), &info);
}
