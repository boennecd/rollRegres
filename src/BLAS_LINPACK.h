/* added to avoid warnings with conflicting definitions in armadillo and R
 * headers. */

void dtrsm(
    const char *side, const char *uplo, const char *transa, const char *diag,
    const int *m, const int *n, const double *alpha, const double *a,
    const int *lda, double *b, const int *ldb);

void factor_first_window(double *x, int rows, int cols, double *qty,
                         double *qraux);

int invert_crossproduct(double *factor, int size);
