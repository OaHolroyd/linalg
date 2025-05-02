/**
 * For the non-cyclic solve this is just the same algorithm as in the
 * pentadiagonal case but with the outer diagonals set to zero.
 *
 * For the cyclic solve we use the Thomas algorithm variant, which makes use of
 * the Sherman-Morrison formula. See
 *   https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm#Variants
 */

#include "tri_solve.h"

#include <string.h>

void tri_lu_factorise(const double *l, double *d, double *u, const int n) {
  /*
   * this is the same as the general LU factorisation without pivoting, but
   * since the entries of A are zero apart from the central diagonals we can
   * unroll the inner loop explicitly.
   */

  // first row is a special case
  u[0] /= d[0];

  // central rows are the same
  for (int i = 1; i < n - 1; i++) {
    d[i] -= l[i] * u[i - 1];
    u[i] /= d[i];
  }

  // last row is a special case
  d[n - 1] -= l[n - 1] * u[n - 2];
}

void tri_lu_solve(
    const double *l, const double *d, const double *u, double *f, const int n
) {
  // solve Ly = f via forward substitution
  double *y = f;
  y[0] /= d[0];
  for (int i = 1; i < n; i++) {
    y[i] = (y[i] - l[i] * y[i - 1]) / d[i];
  }

  // solve Ux = y via backward substitution
  double *x = y;
  for (int i = n - 2; i >= 0; i--) {
    x[i] -= u[i] * x[i + 1];
  }
}

void tri_solve(const double *l, double *d, double *u, double *f, const int n) {
  tri_lu_factorise(l, d, u, n);
  tri_lu_solve(l, d, u, f, n);
}

void cyclic_tri_lu_factorise(
    const double *l, double *d, double *u, double *q, const int n
) {
  // perturb A to get B (also ignoring the periodic entries)
  const double gamma = -d[0];
  d[0] -= gamma;
  d[n - 1] -= u[n - 1] * l[0] / gamma;

  // set g = [gamma, 0, 0, ..., 0, u[n-1]]
  memset(q, 0, n * sizeof(double));
  q[0] = gamma;
  q[n - 1] = u[n - 1];

  // solve q = B \ g, also storing the LU factorisation of B for later
  tri_solve(l, d, u, q, n);
}

void cyclic_tri_lu_solve(
    const double *l, const double *d, const double *u, const double *q,
    double *f, const int n
) {
  // solve y = B \ f
  tri_lu_solve(l, d, u, f, n);

  // compute v[n-1] = l[0] / gamma
  const double gamma = -0.5 * d[0];
  const double vn_1 = l[0] / gamma;

  // compute v·y / (1 + v·q)
  // v·z = z[0] + v[n-1] + z[n-1]
  const double scale =
      (f[0] + vn_1 * f[n - 1]) / (1.0 + q[0] + vn_1 * q[n - 1]);

  // then x = y - q * scale
  for (int i = 0; i < n; i++) {
    f[i] -= q[i] * scale;
  }
}

void cyclic_tri_solve(
    const double *l, double *d, double *u, double *q, double *f, const int n
) {
  cyclic_tri_lu_factorise(l, d, u, q, n);
  cyclic_tri_lu_solve(l, d, u, q, f, n);
}
