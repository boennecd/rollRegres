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
