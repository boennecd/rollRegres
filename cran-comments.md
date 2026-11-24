## Resubmission

rollRegres 0.1.4 was archived on 2023-09-18 for incorrect use of headers.
Version 0.1.5 includes R_ext/BLAS.h directly for dtrsm and the Fortran
character-length arguments. It removes the explicit C++11 specification.
The first window uses the original unpivoted QR routines declared in R_ext/Linpack.h.

## Local check

Ubuntu 24.04.5 LTS, x86_64, R 4.6.2, GCC and GNU Fortran 13.3.0.

R CMD check --as-cran reported zero errors, zero warnings, and two notes.

* Incoming feasibility identifies the package as a new submission and lists
  the archive record.
* The Ubuntu build of R includes -mno-omit-leaf-frame-pointer in its compiler
  flags. The package Makevars does not set this flag.
