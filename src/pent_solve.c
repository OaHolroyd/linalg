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
    const double *l2, double *l1, double *d0, double *u1, double *u2,
    const int n
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
    const double *u2, double *f, const int n
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
    const double *l2, double *l1, double *d0, double *u1, double *u2, double *f,
    const int n
) {
  pent_lu_factorise(l2, l1, d0, u1, u2, n);
  pent_lu_solve(l2, l1, d0, u1, u2, f, n);
}

void cyclic_pent_lu_factorise(
    const double *l2, double *l1, double *d0, double *u1, double *u2,
    double *k0, double *k1, const int n
) {
  // set K = [k0 | k1]
  k0[0] = l2[0];
  for (int i = 1; i < n - 4; i++) {
    k0[i] = 0.0;
  }
  k0[n - 4] = u2[n - 4];
  k0[n - 3] = u1[n - 3];

  k1[0] = l1[0];
  k1[1] = l2[1];
  for (int i = 2; i < n - 3; i++) {
    k1[i] = 0.0;
  }
  k1[n - 3] = u2[n - 3];

  // compute the LU factorisation of E and solve E \ K
  pent_lu_factorise(l2, l1, d0, u1, u2, n - 2);
  pent_lu_solve(l2, l1, d0, u1, u2, k0, n - 2);
  pent_lu_solve(l2, l1, d0, u1, u2, k1, n - 2);

  // compute the 2x2 matrix C - H E^-1 K and store it at the end of the
  // diagonals which have been copied into K.
  d0[n - 2] -=
      u2[n - 2] * k0[0] + l2[n - 2] * k0[n - 4] + l1[n - 2] * k0[n - 3];
  u1[n - 2] -=
      u2[n - 2] * k1[0] + l2[n - 2] * k1[n - 4] + l1[n - 2] * k1[n - 3];
  l1[n - 1] -= u1[n - 1] * k0[0] + u2[n - 1] * k0[1] + l2[n - 1] * k0[n - 3];
  d0[n - 1] -= u1[n - 1] * k1[0] + u2[n - 1] * k1[1] + l2[n - 1] * k1[n - 3];
}

void cyclic_pent_lu_solve(
    const double *l2, const double *l1, const double *l0, const double *u1,
    const double *u2, const double *k0, const double *k1, double *f, const int n
) {
  /*
   * The first step is to solve for the final two elements of the solution:
   *   x[-2:] = (C - H E^-1 K) \ (f[-2:] - H E^-1 f[:-2])
   */

  // solve E \ f[:-2]
  pent_lu_solve(l2, l1, l0, u1, u2, f, n - 2);

  // compute rhs vector for the 2 by 2 subsystem
  f[n - 2] -= u2[n - 2] * f[0] + l2[n - 2] * f[n - 4] + l1[n - 2] * f[n - 3];
  f[n - 1] -= u1[n - 1] * f[0] + u2[n - 1] * f[1] + l2[n - 1] * f[n - 3];

  // solve for the final two elements of the solution vector
  const double det = l0[n - 2] * l0[n - 1] - u1[n - 2] * l1[n - 1];
  const double tmp = (l0[n - 1] * f[n - 2] - u1[n - 2] * f[n - 1]) / det;
  f[n - 1] = (l0[n - 2] * f[n - 1] - l1[n - 1] * f[n - 2]) / det;
  f[n - 2] = tmp;

  /*
   * The complete solution is then found by computing
   *   x[:-2] = E \ f[:-2] - (E \ K) x[-2:]
   */

  // since we've already solved E \ f[:-2] and computed (E \ K) this is
  // relatively simple.
  for (int i = 0; i < n - 2; i++) {
    f[i] -= k0[i] * f[n - 2] + k1[i] * f[n - 1];
  }
}

void cyclic_pent_solve(
    const double *l2, double *l1, double *d0, double *u1, double *u2,
    double *k0, double *k1, double *f, const int n
) {
  cyclic_pent_lu_factorise(l2, l1, d0, u1, u2, k0, k1, n);
  cyclic_pent_lu_solve(l2, l1, d0, u1, u2, k0, k1, f, n);
}
