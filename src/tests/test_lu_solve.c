#include "testing.h"

#include "src/alloc.h"
#include "src/lu_solve.h"

#include <src/io.h>
#include <stdlib.h>

int main(void) {
  START_TEST("lu solve");

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

  /* check LU factorisation solve with multiple RHSs */
  SUBTEST("LU solve multi") {
    const int n = 5;
    const int m = 3;
    double **A = malloc_d2d(n, n);
    double **AA = malloc_d2d(n, n);
    double **F = malloc_d2d(n, m);
    double **FF = malloc_d2d(n, m);
    int *piv = malloc(n * sizeof(int));

    // fill the matrix and rhs with random values
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i][j] = (double)(rand() % 1000 - 500) / 100.0;
        AA[i][j] = A[i][j]; // copy the original matrix
      }

      for (int j = 0; j < m; j++) {
        F[i][j] = (double)(rand() % 1000 - 500) / 100.0;
        FF[i][j] = F[i][j]; // copy the original rhs
      }
    }

    lu_solve_multi(A[0], F[0], piv, n, m);

    // check that AX = F
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        // compute the ijth entry of AX
        double AXij = 0.0;
        for (int k = 0; k < n; k++) {
          AXij += AA[i][k] * F[k][j];
        }
        REQUIRE_CLOSE(AXij, FF[i][j], 1e-10);
      }
    }

    // perform another solve reusing the LU factorisation
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        F[i][j] = (double)(rand() % 1000 - 500) / 100.0;
        FF[i][j] = F[i][j]; // copy the original rhs
      }
    }
    lu_solve_factorised_multi(A[0], piv, F[0], n, m);

    // check that AX = F
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        // compute the ijth entry of AX
        double AXij = 0.0;
        for (int k = 0; k < n; k++) {
          AXij += AA[i][k] * F[k][j];
        }
        REQUIRE_CLOSE(AXij, FF[i][j], 1e-10);
      }
    }

    free_2d(A);
    free_2d(AA);
    free_2d(F);
    free_2d(FF);
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

  /* check LU factorisation without pivoting */
  SUBTEST("LU factorisation no pivot") {
    const int n = 5;
    double **A = malloc_d2d(n, n);
    double **AA = malloc_d2d(n, n);

    // fill the matrix with random values
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i][j] = (double)(rand() % 1000 - 500) / 100.0;
        AA[i][j] = A[i][j]; // copy the original matrix
      }
    }

    lu_factorise(A[0], NULL, n);

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
        REQUIRE_CLOSE(LU[i][j], AA[i][j], 1e-10);
      }
    }

    free_2d(A);
    free_2d(LU);
    free_2d(AA);
  }

  /* check LU factorisation solve without pivoting */
  SUBTEST("LU solve no pivot") {
    const int n = 5;
    double **A = malloc_d2d(n, n);
    double **AA = malloc_d2d(n, n);
    double *f = malloc(n * sizeof(double));
    double *ff = malloc(n * sizeof(double));

    // fill the matrix and rhs with random values
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i][j] = (double)(rand() % 1000 - 500) / 100.0;
        AA[i][j] = A[i][j]; // copy the original matrix
      }
      f[i] = (double)(rand() % 1000 - 500) / 100.0;
      ff[i] = f[i]; // copy the original rhs
    }

    lu_solve(A[0], f, NULL, n);

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
    lu_solve_factorised(A[0], NULL, f, n);

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
  }

  /* check the LU factorisation fails if A is singular without pivoting */
  SUBTEST("LU factorisation singular no pivot") {
    const int n = 5;
    double **A = malloc_d2d(n, n);

    // fill the matrix with random values including a zero row
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i][j] = (i == 0) ? 0.0 : (double)(rand() % 1000 - 500) / 100.0;
      }
    }

    // no pivoting means this should fail immediately
    int err = lu_factorise(A[0], NULL, n);
    REQUIRE(err == 1);

    // fill the matrix with random values including a zero column
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i][j] = (j == 0) ? 0.0 : (double)(rand() % 1000 - 500) / 100.0;
      }
    }

    // as with pivoting, should also fail immediately
    err = lu_factorise(A[0], NULL, n);
    REQUIRE(err == 1);

    free_2d(A);
  }

  /* check the LU factorisation fails if A requires pivoting */
  SUBTEST("LU factorisation requires pivot") {
    const int n = 5;
    double **A = malloc_d2d(n, n);
    double **AA = malloc_d2d(n, n);
    int *piv = malloc(n * sizeof(int));

    // fill the matrix with random values with a zero on the diagonal
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i][j] = (double)(rand() % 1000 - 500) / 100.0;
        AA[i][j] = A[i][j]; // copy the original matrix
      }
    }
    A[0][0] = 0.0;
    AA[0][0] = 0.0;

    // cannot pivot, so a zero on the diagonal is fatal
    int err = lu_factorise(A[0], NULL, n);
    REQUIRE(err == 1);

    // we can get around this with pivoting
    err = lu_factorise(AA[0], piv, n);
    REQUIRE(err == 0);

    free_2d(A);
    free_2d(AA);
    free(piv);
  }

  END_TEST();
}
