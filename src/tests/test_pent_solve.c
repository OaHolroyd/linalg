#include "testing.h"

#include <stdlib.h>
#include <string.h>

#include "src/alloc.h"
#include "src/pent_solve.h"

/**
 * Set the elements of the full, pentadiagonal, matrix A from its diagonals.
 */
static void pent_to_full(
    double **A, double *l2, double *l1, double *d0, double *u1, double *u2,
    int n
) {
  // set all entries to zero
  memset(A[0], 0, n * n * sizeof(double));

  // first row ignores both lower diagonals
  A[0][0] = d0[0];
  A[0][1] = u1[0];
  A[0][2] = u2[0];

  // second row ignores the second lower diagonal
  A[1][0] = l1[1];
  A[1][1] = d0[1];
  A[1][2] = u1[1];
  A[1][3] = u2[1];

  // central rows have all diagonals
  for (int i = 2; i < n - 2; i++) {
    A[i][i - 2] = l2[i];
    A[i][i - 1] = l1[i];
    A[i][i] = d0[i];
    A[i][i + 1] = u1[i];
    A[i][i + 2] = u2[i];
  }

  // penultimate row ignores the first upper diagonal
  A[n - 2][n - 4] = l2[n - 2];
  A[n - 2][n - 3] = l1[n - 2];
  A[n - 2][n - 2] = d0[n - 2];
  A[n - 2][n - 1] = u1[n - 2];

  // last row ignores both upper diagonals
  A[n - 1][n - 3] = l2[n - 1];
  A[n - 1][n - 2] = l1[n - 1];
  A[n - 1][n - 1] = d0[n - 1];
}

/**
 * Set the elements of the full, cyclic, pentadiagonal, matrix A from its
 * diagonals.
 */
static void cyclic_pent_to_full(
    double **A, double *l2, double *l1, double *d0, double *u1, double *u2,
    int n
) {
  // Most elements are as for the non-cyclic case, only the corners need setting
  // separately
  pent_to_full(A, l2, l1, d0, u1, u2, n);

  // first row wraps both lower diagonals
  A[0][n - 2] = l2[0];
  A[0][n - 1] = l1[0];

  // second row wraps the second lower diagonal
  A[1][n - 1] = l2[1];

  // penultimate row wraps the first upper diagonal
  A[n - 2][0] = u2[n - 2];

  // last row wraps both upper diagonals
  A[n - 1][0] = u1[n - 1];
  A[n - 1][1] = u2[n - 1];
}

int main(void) {
  START_TEST("pent solve");

  /* check LU factorisation */
  SUBTEST("pent LU factorisation") {
    const int n = 5;
    double *l2 = malloc(n * sizeof(double));
    double *l1 = malloc(n * sizeof(double));
    double *d0 = malloc(n * sizeof(double));
    double *u1 = malloc(n * sizeof(double));
    double *u2 = malloc(n * sizeof(double));
    double **A = calloc_d2d(n, n);

    // fill the matrix with random values
    for (int i = 0; i < n; i++) {
      l2[i] = (double)(rand() % 1000 - 500) / 100.0;
      l1[i] = (double)(rand() % 1000 - 500) / 100.0;
      d0[i] = (double)(rand() % 1000 - 500) / 100.0;
      u1[i] = (double)(rand() % 1000 - 500) / 100.0;
      u2[i] = (double)(rand() % 1000 - 500) / 100.0;

      // ensure that the matrix is diagonally dominant
      const double d0_mag =
          fabs(d0[i]) + fabs(l1[i]) + fabs(u1[i]) + fabs(l2[i]) + fabs(u2[i]);
      d0[i] = 1.1 * d0_mag * d0[i] / fabs(d0[i]);
    }

    // copy the matrix to a 2D array for checking
    pent_to_full(A, l2, l1, d0, u1, u2, n);

    // factorise and copy to full for checking
    pent_lu_factorise(l2, l1, d0, u1, u2, n);
    double **LU = malloc_d2d(n, n);
    pent_to_full(LU, l2, l1, d0, u1, u2, n);

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

    free(l2);
    free(l1);
    free(d0);
    free(u1);
    free(u2);
    free_2d(A);
    free_2d(LU);
  }

  /* check LU factorisation solve */
  SUBTEST("pent LU solve") {
    const int n = 7;
    double *l2 = malloc(n * sizeof(double));
    double *l1 = malloc(n * sizeof(double));
    double *d0 = malloc(n * sizeof(double));
    double *u1 = malloc(n * sizeof(double));
    double *u2 = malloc(n * sizeof(double));
    double **A = calloc_d2d(n, n);
    double *f = malloc(n * sizeof(double));
    double *ff = malloc(n * sizeof(double));

    // fill the matrix and rhs with random values
    for (int i = 0; i < n; i++) {
      l2[i] = (double)(rand() % 1000 - 500) / 100.0;
      l1[i] = (double)(rand() % 1000 - 500) / 100.0;
      d0[i] = (double)(rand() % 1000 - 500) / 100.0;
      u1[i] = (double)(rand() % 1000 - 500) / 100.0;
      u2[i] = (double)(rand() % 1000 - 500) / 100.0;

      // ensure that the matrix is diagonally dominant
      const double d0_mag =
          fabs(d0[i]) + fabs(l1[i]) + fabs(u1[i]) + fabs(l2[i]) + fabs(u2[i]);
      d0[i] = 1.1 * d0_mag * d0[i] / fabs(d0[i]);

      f[i] = (double)(rand() % 1000 - 500) / 100.0;
      ff[i] = f[i]; // copy the original rhs
    }

    // copy the matrix to a 2D array for checking
    pent_to_full(A, l2, l1, d0, u1, u2, n);

    pent_solve(l2, l1, d0, u1, u2, f, n);

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
    pent_lu_solve(l2, l1, d0, u1, u2, f, n);

    // check that Ax = f
    for (int i = 0; i < n; i++) {
      // compute the ith entry of Ax
      double Axi = 0.0;
      for (int j = 0; j < n; j++) {
        Axi += A[i][j] * f[j];
      }
      REQUIRE_CLOSE(Axi, ff[i], 1e-10);
    }

    free(l2);
    free(l1);
    free(d0);
    free(u1);
    free(u2);
    free_2d(A);
    free(f);
    free(ff);
  }

  /* check LU factorisation solve */
  SUBTEST("cyclic pent LU solve") {
    const int n = 7;
    double *l2 = malloc(n * sizeof(double));
    double *l1 = malloc(n * sizeof(double));
    double *d0 = malloc(n * sizeof(double));
    double *u1 = malloc(n * sizeof(double));
    double *u2 = malloc(n * sizeof(double));
    double *k0 = malloc(n * sizeof(double));
    double *k1 = malloc(n * sizeof(double));
    double **A = calloc_d2d(n, n);
    double *f = malloc(n * sizeof(double));
    double *ff = malloc(n * sizeof(double));

    // fill the matrix and rhs with random values
    for (int i = 0; i < n; i++) {
      l2[i] = (double)(rand() % 1000 - 500) / 100.0;
      l1[i] = (double)(rand() % 1000 - 500) / 100.0;
      d0[i] = (double)(rand() % 1000 - 500) / 100.0;
      u1[i] = (double)(rand() % 1000 - 500) / 100.0;
      u2[i] = (double)(rand() % 1000 - 500) / 100.0;

      // ensure that the matrix is diagonally dominant
      const double d0_mag =
          fabs(d0[i]) + fabs(l1[i]) + fabs(u1[i]) + fabs(l2[i]) + fabs(u2[i]);
      d0[i] = 1.1 * d0_mag * d0[i] / fabs(d0[i]);

      f[i] = (double)(rand() % 1000 - 500) / 100.0;
      ff[i] = f[i]; // copy the original rhs
    }

    // copy the matrix to a 2D array for checking
    cyclic_pent_to_full(A, l2, l1, d0, u1, u2, n);

    cyclic_pent_solve(l2, l1, d0, u1, u2, k0, k1, f, n);

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
    cyclic_pent_lu_solve(l2, l1, d0, u1, u2, k0, k1, f, n);

    // check that Ax = f
    for (int i = 0; i < n; i++) {
      // compute the ith entry of Ax
      double Axi = 0.0;
      for (int j = 0; j < n; j++) {
        Axi += A[i][j] * f[j];
      }
      REQUIRE_CLOSE(Axi, ff[i], 1e-10);
    }

    free(l2);
    free(l1);
    free(d0);
    free(u1);
    free(u2);
    free(k0);
    free(k1);
    free_2d(A);
    free(f);
    free(ff);
  }

  END_TEST();
}
