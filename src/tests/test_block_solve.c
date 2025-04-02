#include "testing.h"

#include "src/alloc.h"
#include "src/block_solve.h"
#include "src/lu_solve.h"

#include <src/io.h>
#include <stdlib.h>

int main(void) {
  START_TEST("block solve");

  /* check LU factorisation */
  SUBTEST("block solve") {
    const int n = 5;
    const int m = 5;
    const int nm = n + m;
    double **R = malloc_d2d(nm, nm);
    int *piv = malloc(nm * sizeof(int));
    double *f = malloc(nm * sizeof(double));
    double **A = malloc_d2d(n, n);
    double **B = malloc_d2d(n, m);
    double **C = malloc_d2d(m, n);
    double **D = malloc_d2d(m, m);
    double *a = malloc(n * sizeof(double));
    double *b = malloc(m * sizeof(double));

    // fill the matrix with random values
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        R[i][j] = (double)(rand() % 1000 - 500) / 100.0;
        A[i][j] = R[i][j];
      }
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        R[i][n + j] = (double)(rand() % 1000 - 500) / 100.0;
        B[i][j] = R[i][n + j];
      }
    }
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        R[n + i][j] = (double)(rand() % 1000 - 500) / 100.0;
        C[i][j] = R[n + i][j];
      }
    }
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
        R[n + i][n + j] = (double)(rand() % 1000 - 500) / 100.0;
        D[i][j] = R[n + i][n + j];
      }
    }

    // fill the rhs with random values
    for (int i = 0; i < n; i++) {
      f[i] = (double)(rand() % 1000 - 500) / 100.0;
      a[i] = f[i];
    }
    for (int i = 0; i < m; i++) {
      f[n + i] = (double)(rand() % 1000 - 500) / 100.0;
      b[i] = f[n + i];
    }

    int err = lu_solve(R[0], f, piv, nm);
    REQUIRE_BARRIER(err == 0);

    size_t work_size = n * m + ((n > m) ? n : m);
    double *work = malloc(work_size * sizeof(double));
    int *pivn = piv;
    int *pivm = piv + n;
    err = block_solve(A[0], B[0], C[0], D[0], a, b, pivn, pivm, work, n, m);
    REQUIRE_BARRIER(err == 0);

    // check that the block solve and the normal solve match
    for (int i = 0; i < n; i++) {
      REQUIRE_CLOSE(f[i], a[i], 1e-10);
    }
    for (int i = 0; i < m; i++) {
      REQUIRE_CLOSE(f[n + i], b[i], 1e-10);
    }

    free_2d(R);
    free_2d(A);
    free_2d(B);
    free_2d(C);
    free_2d(D);
    free(piv);
    free(f);
    free(a);
    free(b);
    free(work);
  }

  END_TEST();
}
