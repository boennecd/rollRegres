#include "BLAS_LINPACK.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include "R_ext/Linpack.h"

void dtrsm(
    const char *side, const char *uplo, const char *transa, const char *diag,
    const int *m, const int *n, const double *alpha, const double *a,
    const int *lda, double *b, const int *ldb){
  return F77_CALL(dtrsm)(
      side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb
      FCONE FCONE FCONE FCONE);
}

void factor_first_window(double *x, int rows, int cols, double *qty,
                         double *qraux){
  int no_pivot = 0, unused_pivot = 0;
  double unused = 0.;
  F77_CALL(dqrdc)(x, &rows, &rows, &cols, qraux, &unused_pivot, &unused,
                  &no_pivot);

  // Job 1000 only computes Q^T y. LINPACK permits y and qty to coincide.
  int job = 1000, info = 0;
  F77_CALL(dqrsl)(x, &rows, &rows, &cols, qraux, qty, &unused, qty,
                  &unused, &unused, &unused, &job, &info);
}

int invert_crossproduct(double *factor, int size){
  const char upper = 'U';
  int info = 0;
  F77_CALL(dpotri)(&upper, &size, factor, &size, &info FCONE);
  return info;
}
