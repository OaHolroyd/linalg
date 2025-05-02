/**
 * The standard pentadiagonal solver is described in 'Numerical Recipes in C',
 * also described in section D of 'cuPentBatch – A batched pentadiagonal solver
 * for NVIDIA GPUs': https://arxiv.org/pdf/1807.07382.
 *
 * The periodic solver is described in section C of the same paper, but is
 * originally described in 'Pent: A periodic pentadiagonal systems solver'.
 */

#include "pent_solve.h"

void pent_lu_factorise(
    const double *l2, double *l1, double *d0, double *u1, double *u2, int n
) {
  /*
   * this is the same as the general LU factorisation without pivoting, but
   * since the entries of A are zero apart from the central diagonals we can
   * unroll the inner loop explicitly.
   */

  // first row is a special case
  u1[0] /= d0[0];
  u2[0] /= d0[0];

  // second row is a special case
  d0[1] -= l1[1] * u1[0];
  u1[1] = (u1[1] - l1[1] * u2[0]) / d0[1];
  u2[1] /= d0[1];

  // central rows are the same
  for (int i = 2; i < n - 2; i++) {
    l1[i] -= l2[i] * u1[i - 2];
    d0[i] -= l2[i] * u2[i - 2] + l1[i] * u1[i - 1];
    u1[i] = (u1[i] - l1[i] * u2[i - 1]) / d0[i];
    u2[i] /= d0[i];
  }

  // penultimate row is a special case
  l1[n - 2] -= l2[n - 2] * u1[n - 4];
  d0[n - 2] -= l2[n - 2] * u2[n - 4] + l1[n - 2] * u1[n - 3];
  u1[n - 2] = (u1[n - 2] - l1[n - 2] * u2[n - 3]) / d0[n - 2];

  // last row is a special case
  l1[n - 1] -= l2[n - 1] * u1[n - 3];
  d0[n - 1] -= l2[n - 1] * u2[n - 3] + l1[n - 1] * u1[n - 2];
}

void pent_lu_solve(
    const double *l2, const double *l1, const double *l0, const double *u1,
    const double *u2, double *f, int n
) {
  // solve Ly = f via forward substitution
  double *y = f;
  y[0] /= l0[0];
  y[1] = (y[1] - l1[1] * y[0]) / l0[1];
  for (int i = 2; i < n; i++) {
    y[i] = (y[i] - l1[i] * y[i - 1] - l2[i] * y[i - 2]) / l0[i];
  }

  // solve Ux = y via backward substitution
  double *x = y;
  x[n - 2] -= u1[n - 2] * x[n - 1];
  for (int i = n - 3; i >= 0; i--) {
    x[i] -= u1[i] * x[i + 1] + u2[i] * x[i + 2];
  }
}

void pent_solve(
    double *l2, double *l1, double *d0, double *u1, double *u2, double *f, int n
) {
  pent_lu_factorise(l2, l1, d0, u1, u2, n);
  pent_lu_solve(l2, l1, d0, u1, u2, f, n);
}
