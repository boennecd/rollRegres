## New version submission

Version 0.1.6 adds optional coefficient covariance matrices and standard
errors for rolling and expanding regressions.

CRAN archived rollRegres 0.1.4 on 2023-09-18 for incorrect use of headers.
The package now includes `R_ext/BLAS.h` directly for `dtrsm` and its Fortran
character-length arguments.

## Local check

Debian GNU/Linux 13, x86_64, R 4.5.3, GCC and GNU Fortran 16.2.0.

`R CMD check --as-cran rollRegres_0.1.6.tar.gz` reported zero errors, zero
warnings, and one note. The incoming check lists the archive record.

The local check removed the Conda toolchain's target-specific `-march` flag.
It skipped the online clock check because the time service was unavailable.
