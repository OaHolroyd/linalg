/**
 * This block solve is described in "Inverses of 2x2 Block Matrices" by Lu and
 * Shiou 2000 (https://doi.org/10.1016/S0898-1221(01)00278-4), eqn 2.2.
 */

#include "block_solve.h"

#include <string.h>

#include "lu_solve.h"

int block_solve(
    double *A, const double *B, const double *C, double *D, double *a,
    double *b, int *pivn, int *pivm, double *work, const int n, const int m
) {
  // solve a = A \ a
  int err = lu_solve(A, a, pivn, n);
  if (err != 0) {
    return err;
  }

  // store a copy of B
  double *AB = work;
  memcpy(AB, B, m * n * sizeof(double));

  // compute S = D - C A \ B
  lu_solve_factorised_multi(A, pivn, AB, n, m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      for (int k = 0; k < n; k++) {
        D[i * m + j] -= C[i * n + k] * AB[k * m + j];
      }
    }
  }

  // compute z = C a
  double *z = work + m * n; // here we're using m entries, later we will use n
  for (int i = 0; i < m; i++) {
    z[i] = 0.0;
    for (int j = 0; j < n; j++) {
      z[i] += C[i * n + j] * a[j];
    }
  }

  // solve z = S \ z, b = S \ b
  err = lu_factorise(D, pivm, m);
  if (err != 0) {
    return err;
  }
  lu_solve_factorised(D, pivm, z, m);
  lu_solve_factorised(D, pivm, b, m);

  // compute b = b - z
  for (int i = 0; i < m; i++) {
    b[i] -= z[i];
  }

  // z = B b
  // now using m entries of z
  for (int i = 0; i < n; i++) {
    z[i] = 0.0;
    for (int j = 0; j < m; j++) {
      z[i] += B[i * m + j] * b[j];
    }
  }

  // z = A \ z
  lu_solve_factorised(A, pivn, z, n);

  // a = a - z
  for (int i = 0; i < n; i++) {
    a[i] -= z[i];
  }

  return 0;
}

int block_solve_simplified(
    const double *B, const double *C, double *S, double *a, double *b,
    int *pivm, double *work, const int n, const int m
) {
  // compute z = C a
  double *z = work; // here we're using m entries, later we will use n
  for (int i = 0; i < m; i++) {
    z[i] = 0.0;
    for (int j = 0; j < n; j++) {
      z[i] += C[i * n + j] * a[j];
    }
  }

  // solve z = S \ z, b = S \ b
  const int err = lu_factorise(S, pivm, m);
  if (err != 0) {
    return err;
  }
  lu_solve_factorised(S, pivm, z, m);
  lu_solve_factorised(S, pivm, b, m);

  // compute b = b - z
  for (int i = 0; i < m; i++) {
    b[i] -= z[i];
  }

  // z = B b
  // now using m entries of z
  for (int i = 0; i < n; i++) {
    z[i] = 0.0;
    for (int j = 0; j < m; j++) {
      z[i] += B[i * m + j] * b[j];
    }
  }

  // a = a - z
  for (int i = 0; i < n; i++) {
    a[i] -= z[i];
  }

  return 0;
}
