#include "BLAS_LINPACK.h"
#include "R_ext/Linpack.h"

void dtrsm(
    const char *side, const char *uplo, const char *transa, const char *diag,
    const int *m, const int *n, const double *alpha, const double *a,
    const int *lda, double *b, const int *ldb){
  return F77_CALL(dtrsm)(
      side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb
      FCONE FCONE FCONE FCONE);
}
