#include "lu-solve.h"

#include <math.h>

#define LU_TOL (1e-10)

int lu_factorise(double *A, int *piv, const int n) {
  // start with a unit pivot matrix
  for (int i = 0; i < n; i++) {
    piv[i] = i;
  }

  for (int i = 0; i < n; i++) {
    double maxA = 0.0;
    int maxi = i;

    // find the largest entry in the column not above the diagonal
    for (int j = i; j < n; j++) {
      if (fabs(A[j * n + i]) > maxA) {
        maxA = fabs(A[j * n + i]);
        maxi = j;
      }
    }

    // if the largest entry is too small, the matrix is singular
    if (maxA < LU_TOL) {
      return i + 1; // return the row of the first zero pivot
    }

    // pivot if necessary
    if (maxi != i) {
      // swap the rows in the pivot array
      const int tmp = piv[i];
      piv[i] = piv[maxi];
      piv[maxi] = tmp;

      // swap the rows in the matrix
      for (int j = 0; j < n; j++) {
        const double tmpA = A[i * n + j];
        A[i * n + j] = A[maxi * n + j];
        A[maxi * n + j] = tmpA;
      }
    }

    for (int j = i + 1; j < n; j++) {
      // divide the pivot row by the pivot element
      A[j * n + i] /= A[i * n + i];

      // subtract the pivot row from the current row (Gaussian elimination)
      for (int k = i + 1; k < n; k++) {
        A[j * n + k] -= A[j * n + i] * A[i * n + k];
      }
    }
  }

  return 0;
}

int lu_solve_factorised(
    const double *LU, const int *piv, double *f, const int n
) {
  // solve Ly = Pf by forward substitution
  for (int i = 0; i < n; i++) {
    // pivot the right-hand side
    if (i < piv[i]) {
      const double tmp = f[i];
      f[i] = f[piv[i]];
      f[piv[i]] = tmp;
    }

    for (int k = 0; k < i; k++) {
      f[i] -= LU[i * n + k] * f[k];
    }
  }

  // solve Ux = y by back substitution
  for (int i = n - 1; i >= 0; i--) {
    for (int k = i + 1; k < n; k++) {
      f[i] -= LU[i * n + k] * f[k];
    }
    f[i] /= LU[i * n + i];
  }

  return 0;
}

int lu_solve(double *A, double *f, int *piv, int n) {
  // factorise the matrix
  int err = lu_factorise(A, piv, n);
  if (err != 0) {
    return err; // return the row of the first zero pivot
  }

  // solve the factorised system of equations
  return lu_solve_factorised(A, piv, f, n);
}
