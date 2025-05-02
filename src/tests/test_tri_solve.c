#include "testing.h"

#include <stdlib.h>
#include <string.h>

#include "src/alloc.h"
#include "src/tri_solve.h"

/**
 * Set the elements of the full, tridiagonal, matrix A from its diagonals.
 */
static void tri_to_full(double **A, double *l, double *d, double *u, int n) {
  // set all entries to zero
  memset(A[0], 0, n * n * sizeof(double));

  // first row ignores lower diagonal
  A[0][0] = d[0];
  A[0][1] = u[0];

  // central rows have all diagonals
  for (int i = 1; i < n - 1; i++) {
    A[i][i - 1] = l[i];
    A[i][i] = d[i];
    A[i][i + 1] = u[i];
  }

  // last row ignores upper diagonal
  A[n - 1][n - 2] = l[n - 1];
  A[n - 1][n - 1] = d[n - 1];
}

/**
 * Set the elements of the full, cyclic, tridiagonal, matrix A from its
 * diagonals.
 */
static void
cyclic_tri_to_full(double **A, double *l, double *d, double *u, int n) {
  // Most elements are as for the non-cyclic case, only the corners need setting
  // separately
  tri_to_full(A, l, d, u, n);

  // first row wraps lower diagonal
  A[0][n - 1] = l[0];

  // last row wraps upper diagonal
  A[n - 1][0] = u[n - 1];
}

int main(void) {
  START_TEST("tri solve");

  /* check LU factorisation */
  SUBTEST("tri LU factorisation") {
    const int n = 5;
    double *l = malloc(n * sizeof(double));
    double *d = malloc(n * sizeof(double));
    double *u = malloc(n * sizeof(double));
    double **A = calloc_d2d(n, n);

    // fill the matrix with random values
    for (int i = 0; i < n; i++) {
      l[i] = (double)(rand() % 1000 - 500) / 100.0;
      d[i] = (double)(rand() % 1000 - 500) / 100.0;
      u[i] = (double)(rand() % 1000 - 500) / 100.0;

      // ensure that the matrix is diagonally dominant
      const double d0_mag = fabs(d[i]) + fabs(l[i]) + fabs(u[i]);
      d[i] = 1.1 * d0_mag * d[i] / fabs(d[i]);
    }

    // copy the matrix to a 2D array for checking
    tri_to_full(A, l, d, u, n);

    // factorise and copy to full for checking
    tri_lu_factorise(l, d, u, n);
    double **LU = malloc_d2d(n, n);
    tri_to_full(LU, l, d, u, n);

    // multiply the L and U matrices to check that they are correct
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        // compute the ijth element
        double LUij = 0.0;
        for (int k = 0; k < n; k++) {
          // extract the elements of L and U
          const double Lik = (k > i) ? 0.0 : LU[i][k];
          const double Ukj = (k == j) ? 1.0 : ((j < k) ? 0.0 : LU[k][j]);
          LUij += Lik * Ukj;
        }

        REQUIRE_CLOSE(LUij, A[i][j], 1e-10);
      }
    }

    free(l);
    free(d);
    free(u);
    free_2d(A);
    free_2d(LU);
  }

  /* check LU factorisation solve */
  SUBTEST("tri LU solve") {
    const int n = 7;
    double *l = malloc(n * sizeof(double));
    double *d = malloc(n * sizeof(double));
    double *u = malloc(n * sizeof(double));
    double **A = calloc_d2d(n, n);
    double *f = malloc(n * sizeof(double));
    double *ff = malloc(n * sizeof(double));

    // fill the matrix and rhs with random values
    for (int i = 0; i < n; i++) {
      l[i] = (double)(rand() % 1000 - 500) / 100.0;
      d[i] = (double)(rand() % 1000 - 500) / 100.0;
      u[i] = (double)(rand() % 1000 - 500) / 100.0;

      // ensure that the matrix is diagonally dominant
      const double d_mag = fabs(d[i]) + fabs(l[i]) + fabs(u[i]);
      d[i] = 1.1 * d_mag * d[i] / fabs(d[i]);

      f[i] = (double)(rand() % 1000 - 500) / 100.0;
      ff[i] = f[i]; // copy the original rhs
    }

    // copy the matrix to a 2D array for checking
    tri_to_full(A, l, d, u, n);

    tri_solve(l, d, u, f, n);

    // check that Ax = f
    for (int i = 0; i < n; i++) {
      // compute the ith entry of Ax
      double Axi = 0.0;
      for (int j = 0; j < n; j++) {
        Axi += A[i][j] * f[j];
      }
      REQUIRE_CLOSE(Axi, ff[i], 1e-10);
    }

    // perform another solve reusing the LU factorisation
    for (int i = 0; i < n; i++) {
      f[i] = (double)(rand() % 1000 - 500) / 100.0;
      ff[i] = f[i]; // copy the original rhs
    }
    tri_lu_solve(l, d, u, f, n);

    // check that Ax = f
    for (int i = 0; i < n; i++) {
      // compute the ith entry of Ax
      double Axi = 0.0;
      for (int j = 0; j < n; j++) {
        Axi += A[i][j] * f[j];
      }
      REQUIRE_CLOSE(Axi, ff[i], 1e-10);
    }

    free(l);
    free(d);
    free(u);
    free_2d(A);
    free(f);
    free(ff);
  }

  /* check LU factorisation solve */
  SUBTEST("cyclic tri LU solve") {
    const int n = 7;
    double *l = malloc(n * sizeof(double));
    double *d = malloc(n * sizeof(double));
    double *u = malloc(n * sizeof(double));
    double *q = malloc(n * sizeof(double));
    double **A = calloc_d2d(n, n);
    double *f = malloc(n * sizeof(double));
    double *ff = malloc(n * sizeof(double));

    // fill the matrix and rhs with random values
    for (int i = 0; i < n; i++) {
      l[i] = (double)(rand() % 1000 - 500) / 100.0;
      d[i] = (double)(rand() % 1000 - 500) / 100.0;
      u[i] = (double)(rand() % 1000 - 500) / 100.0;

      // ensure that the matrix is diagonally dominant
      const double d_mag = fabs(d[i]) + fabs(l[i]) + fabs(u[i]);
      d[i] = 1.1 * d_mag * d[i] / fabs(d[i]);

      f[i] = (double)(rand() % 1000 - 500) / 100.0;
      ff[i] = f[i]; // copy the original rhs
    }

    // copy the matrix to a 2D array for checking
    cyclic_tri_to_full(A, l, d, u, n);

    cyclic_tri_solve(l, d, u, q, f, n);

    // check that Ax = f
    for (int i = 0; i < n; i++) {
      // compute the ith entry of Ax
      double Axi = 0.0;
      for (int j = 0; j < n; j++) {
        Axi += A[i][j] * f[j];
      }
      REQUIRE_CLOSE(Axi, ff[i], 1e-10);
    }

    // perform another solve reusing the LU factorisation
    for (int i = 0; i < n; i++) {
      f[i] = (double)(rand() % 1000 - 500) / 100.0;
      ff[i] = f[i]; // copy the original rhs
    }
    cyclic_tri_lu_solve(l, d, u, q, f, n);

    // check that Ax = f
    for (int i = 0; i < n; i++) {
      // compute the ith entry of Ax
      double Axi = 0.0;
      for (int j = 0; j < n; j++) {
        Axi += A[i][j] * f[j];
      }
      REQUIRE_CLOSE(Axi, ff[i], 1e-10);
    }

    free(l);
    free(d);
    free(u);
    free(q);
    free_2d(A);
    free(f);
    free(ff);
  }

  END_TEST();
}
