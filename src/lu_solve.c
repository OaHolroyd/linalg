#include "lu_solve.h"

#include <math.h>
#include <stddef.h>

#define LU_TOL (1e-10)

/**
 * Permute the right-hand side vectors according to the permutation array.
 *
 * @param F right-hand side vectors, overwritten with permuted vectors
 * @param piv permutation array of length n, left unchanged
 * @param n number of rows in F
 * @param m number of right-hand side vectors
 */
static void permute_vectors(double *F, int *piv, const int n, const int m) {
  for (int k = 0; k < n; k++) {
    // find the first element that is not in the right place
    int i = k;
    for (; i < n; i++) {
      if (piv[i] >= 0) {
        break;
      }
    }
    if (i == n) {
      break; // all elements are in the right place
    }

    // go through the cycle starting with i until we get back to the start
    const int ii = i; // remember the start of the cycle
    int pi = piv[i];
    piv[i] = -piv[i] - 1; // mark as used
    while (pi != ii) {
      // swap the rows of the right-hand side
      for (int j = 0; j < m; j++) {
        const double tmp = F[i * m + j];
        F[i * m + j] = F[pi * m + j];
        F[pi * m + j] = tmp;
      }

      // move forwards in the cycle
      i = pi;
      pi = piv[pi];
      piv[i] = -piv[i] - 1; // mark as used
    }
  }

  // put piv back to the original form so it can be reused
  for (int i = 0; i < n; i++) {
    piv[i] = -piv[i] - 1;
  }
}

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

void lu_solve_factorised(const double *LU, int *piv, double *f, const int n) {
  // pivot the right-hand side to compute Pf
  if (piv != NULL) {
    permute_vectors(f, piv, n, 1);
  }

  // solve Ly = Pf by forward substitution
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

void lu_solve_factorised_multi(
    const double *LU, int *piv, double *F, int n, int m
) {
  // pivot the right-hand side to compute PF
  if (piv != NULL) {
    permute_vectors(F, piv, n, m);
  }

  // solve LY = PF by forward substitution
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < i; k++) {
      for (int j = 0; j < m; j++) { // apply to entire row
        F[i * m + j] -= LU[i * n + k] * F[k * m + j];
      }
    }
  }

  // solve UX = Y by back substitution
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

int lu_solve(double *A, double *f, int *piv, int n) {
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
