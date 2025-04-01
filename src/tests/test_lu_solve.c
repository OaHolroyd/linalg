#include "testing.h"

#include "src/alloc.h"
#include "src/lu-solve.h"

#include <stdlib.h>

int main(void) {
  START_TEST("alloc");

  /* check LU factorisation */
  SUBTEST("LU factorisation") {
    const int n = 5;
    double **A = malloc_d2d(n, n);
    double **AA = malloc_d2d(n, n);
    int *piv = malloc(n * sizeof(int));

    // fill the matrix with random values
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i][j] = (double)(rand() % 1000 - 500) / 100.0;
        AA[i][j] = A[i][j]; // copy the original matrix
      }
    }

    lu_factorise(A[0], piv, n);

    // multiply the L and U matrices to check that they are correct
    double **LU = malloc_d2d(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        LU[i][j] = 0.0;
        for (int k = 0; k < n; k++) {
          // extract the elements of L and U
          double Lik = (i == k) ? 1.0 : (((i < k) ? 0.0 : A[i][k]));
          double Ukj = (k > j) ? 0.0 : A[k][j];
          LU[i][j] += Lik * Ukj;
        }
      }
    }

    // check that PA = LU
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        REQUIRE_CLOSE(LU[i][j], AA[piv[i]][j], 1e-10);
      }
    }

    free_2d(A);
    free_2d(LU);
    free_2d(AA);
    free(piv);
  }

  /* check LU factorisation solve */
  SUBTEST("LU solve") {
    const int n = 5;
    double **A = malloc_d2d(n, n);
    double **AA = malloc_d2d(n, n);
    double *f = malloc(n * sizeof(double));
    double *ff = malloc(n * sizeof(double));
    int *piv = malloc(n * sizeof(int));

    // fill the matrix and rhs with random values
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i][j] = (double)(rand() % 1000 - 500) / 100.0;
        AA[i][j] = A[i][j]; // copy the original matrix
      }
      f[i] = (double)(rand() % 1000 - 500) / 100.0;
      ff[i] = f[i]; // copy the original rhs
    }

    lu_solve(A[0], f, piv, n);

    // check that Ax = f
    for (int i = 0; i < n; i++) {
      // compute the ith entry of Ax
      double Axi = 0.0;
      for (int j = 0; j < n; j++) {
        Axi += AA[i][j] * f[j];
      }
      REQUIRE_CLOSE(Axi, ff[i], 1e-10);
    }

    // perform another solve reusing the LU factorisation
    for (int i = 0; i < n; i++) {
      f[i] = (double)(rand() % 1000 - 500) / 100.0;
      ff[i] = f[i]; // copy the original rhs
    }
    lu_solve_factorised(A[0], piv, f, n);

    // check that Ax = f
    for (int i = 0; i < n; i++) {
      // compute the ith entry of Ax
      double Axi = 0.0;
      for (int j = 0; j < n; j++) {
        Axi += AA[i][j] * f[j];
      }
      REQUIRE_CLOSE(Axi, ff[i], 1e-10);
    }

    free_2d(A);
    free_2d(AA);
    free(f);
    free(ff);
    free(piv);
  }

  /* check the LU factorisation fails if A is singular */
  SUBTEST("LU factorisation singular") {
    const int n = 5;
    double **A = malloc_d2d(n, n);
    int *piv = malloc(n * sizeof(int));

    // fill the matrix with random values including a zero row
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i][j] = (i == 0) ? 0.0 : (double)(rand() % 1000 - 500) / 100.0;
      }
    }

    // the pivoting should push the zeros down to the final row, resulting
    // in a failure at the last row
    int err = lu_factorise(A[0], piv, n);
    REQUIRE(err == n);

    // fill the matrix with random values including a zero column
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i][j] = (j == 0) ? 0.0 : (double)(rand() % 1000 - 500) / 100.0;
      }
    }

    // pivoting can't do anything with a zero column, so the first row fails
    err = lu_factorise(A[0], piv, n);
    REQUIRE(err == 1);

    free_2d(A);
    free(piv);
  }

  END_TEST();
}
