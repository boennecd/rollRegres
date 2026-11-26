# rollRegres 0.1.6

* Use `do_compute = "vcov"` to compute coefficient covariance matrices for
  rolling and expanding windows. Use `vcov(fit)` to get the covariance array
  or `vcov(fit, se = TRUE)` to get the standard errors.
* Use `coef(fit)` to get the coefficient matrix from either fit interface.
* When the residual norm drops sharply, the fitter recomputes the QR
  decomposition to reduce roundoff from a large residual leaving the window.

# rollRegres 0.1.5

* The package includes the R BLAS header directly. This fixes compilation on
  versions of R that no longer include it through the LINPACK header.
* Rolling regressions now retain valid coefficients when roundoff prevents a
  residual norm downdate. The fitter recomputes the norm when it is needed.
* The formula interface now passes `min_obs` to the fitter. Ungrouped fits
  accept a window width equal to the number of coefficients.
* The fitter reports an error if a fitted window has fewer rows than coefficients.
* The package uses the default C++ standard from R.

# rollRegres 0.1.4
* added missing char array length arguments to Fortran calls.
* fix an issue with the `group` argument when the first group(s) is ignored 
  there is a "hole" in the subsequent window. This previously caused a crash.

# rollRegres 0.1.3
* added ORCID and used `rng = false` in `Rcpp::export` attributes will have 
  reduced the computation time.

# rollRegres 0.1.2
* the c++ code has been re-written. It now does essentially as described in 
  the LINPACK user guide chapter 8 and 9. This is faster and more precise.

# rollRegres 0.1.1
* handle data sets with gaps
* add the option to require a minimum number of observations in a window
* output with groups have changed such that the first window max `grp` value will have 
at least `width - 1L` distance from `min(grp)` in `roll_regres` and 
`roll_regres.fit`
