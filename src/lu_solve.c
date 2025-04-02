#include "lu_solve.h"

#include <math.h>
#include <stddef.h>

#define LU_TOL (1e-10)

int lu_factorise(double *A, int *piv, const int n) {
  if (piv == NULL) {
    return lu_factorise_no_pivoting(A, n);
  }

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

int lu_factorise_no_pivoting(double *A, const int n) {
  for (int i = 0; i < n; i++) {
    // if the diagonal entry is too small, the matrix is singular or requires
    // pivoting to factorise
    if (fabs(A[i * n + i]) < LU_TOL) {
      return i + 1; // return the row of the first zero pivot
    }

    for (int j = i + 1; j < n; j++) {
      // divide the row by the diagonal entry
      A[j * n + i] /= A[i * n + i];

      // subtract the row from the current row (Gaussian elimination)
      for (int k = i + 1; k < n; k++) {
        A[j * n + k] -= A[j * n + i] * A[i * n + k];
      }
    }
  }

  return 0;
}

void lu_solve_factorised(
    const double *LU, const int *piv, double *f, const int n
) {
  if (piv == NULL) {
    lu_solve_factorised_no_pivoting(LU, f, n);
    return;
  }

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
}

void lu_solve_factorised_multi(
    const double *LU, const int *piv, double *F, int n, int m
) {
  // solve Ly = Pf by forward substitution
  for (int i = 0; i < n; i++) {
    // pivot the right-hand side
    if (i < piv[i]) {
      for (int j = 0; j < m; j++) { // apply to entire row
        const double tmp = F[i * m + j];
        F[i * m + j] = F[piv[i] * m + j];
        F[piv[i] * m + j] = tmp;
      }
    }

    for (int k = 0; k < i; k++) {
      for (int j = 0; j < m; j++) { // apply to entire row
        F[i * m + j] -= LU[i * n + k] * F[k * m + j];
      }
    }
  }

  // solve Ux = y by back substitution
  for (int i = n - 1; i >= 0; i--) {
    for (int k = i + 1; k < n; k++) {
      for (int j = 0; j < m; j++) { // apply to entire row
        F[i * m + j] -= LU[i * n + k] * F[k * m + j];
      }
    }
    for (int j = 0; j < m; j++) { // apply to entire row
      F[i * m + j] /= LU[i * n + i];
    }
  }
}

void lu_solve_factorised_no_pivoting(const double *LU, double *f, const int n) {
  // solve Ly = f by forward substitution
  for (int i = 0; i < n; i++) {
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
}

int lu_solve(double *A, double *f, int *piv, int n) {
  if (piv == NULL) {
    return lu_solve_no_pivoting(A, f, n);
  }

  // factorise the matrix
  int err = lu_factorise(A, piv, n);
  if (err != 0) {
    return err; // return the row of the first zero pivot
  }

  // solve the factorised system of equations
  lu_solve_factorised(A, piv, f, n);
  return 0;
}

int lu_solve_multi(double *A, double *F, int *piv, int n, int m) {
  // factorise the matrix
  int err = lu_factorise(A, piv, n);
  if (err != 0) {
    return err; // return the row of the first zero pivot
  }

  // solve the factorised system of equations
  lu_solve_factorised_multi(A, piv, F, n, m);
  return 0;
}

int lu_solve_no_pivoting(double *A, double *f, int n) {
  // factorise the matrix
  int err = lu_factorise_no_pivoting(A, n);
  if (err != 0) {
    return err; // return the row of the first zero pivot
  }

  // solve the factorised system of equations
  lu_solve_factorised_no_pivoting(A, f, n);
  return 0;
}
