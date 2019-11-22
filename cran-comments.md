## Test environments
* Ubuntu 18.04 LTS with gcc 8.3.0
  R version 3.6.1
* Ubuntu 16.04 LTS (on travis-ci)
  R version 3.6.1
* Ubuntu 18.04 LTS with gcc 8.3.0 with --enable-lto
  R devel
* Ubuntu 18.04 LTS with clang 6.0.0 with ASAN and 
  UBSAN checks
  R devel
* win-builder (devel and release)
* `rhub::check_for_cran()` and `rhub::check_on_solaris()`
  
## R CMD check results
There were no ERRORs, WARNINGs, or NOTES.

I have solved the issue with the for misrepresentation of authorship. Now, 
I give credit to the included LINPACK routines and to Madeleine Thompson who
has modified them.

I am sorry for the error in the previous version. I must have been way too
quick when modified the comments in the Fortran files.
