#ifndef TRI_SOLVE_H
#define TRI_SOLVE_H

/**
 * Factorises a tridiagonal, diagonally dominant, square matrix A into A = LU.
 *
 * This is a special case of the LU factorisation, where the matrix A has three
 * diagonals (one lower, one main, one upper) and is diagonally dominant (i.e.
 * the absolute value of the diagonal element is greater than the sum of the
 * absolute values of the other elements in the same row). The resulting matrix
 * L is empty except for the main diagonal and the lower diagonal, and the
 * matrix U is empty except for the main diagonal and the upper diagonal. Since
 * the main diagonal elements of U are all 1, it is not necessary to store them,
 * so the factorisation can be performed in-place. This takes O(n) steps, a
 * significant improvement over the general O(n^3) LU factorisation.
 *
 * @param l lower diagonal, overwritten with lower diagonal of L
 * @param d main diagonal, overwritten with main diagonal of L
 * @param u upper diagonal, overwritten with upper diagonal of U
 * @param n size of the matrix
 */
void tri_lu_factorise(const double *l, double *d, double *u, int n);

/**
 * Given an LU factorisation of a tridiagonal, square, matrix A = LU, solves
 * Ax = f in place.
 *
 * Unlike in the general case (which is O(n^2)), this can be done in O(n) steps.
 *
 * @param l lower diagonal of L
 * @param d main diagonal of L
 * @param u upper diagonal of U
 * @param f right-hand side vector, overwritten with the solution
 * @param n size of the matrix
 */
void tri_lu_solve(
    const double *l, const double *d, const double *u, double *f, int n
);

/**
 * Solves the system Ax = f in place, where A is a tridiagonal.
 *
 * The index into the diagonal represents the row:
 *    d[0]  u[0]     0     0     0     0 ...
 *    l[1]  d[1]  u[1]     0     0     0 ...
 *       0  l[2]  d[2]  u[2]     0     0 ...
 *       0     0  l[3]  d[3]  u[3]     0 ...
 * and so l[0], u[n-1] are not accessed.
 *
 * @param l lower diagonal, overwritten with lower diagonal of L
 * @param d main diagonal, overwritten with main diagonal of L
 * @param u upper diagonal, overwritten with upper diagonal of U
 * @param f right-hand side vector, overwritten with the solution
 * @param n size of the matrix
 */
void tri_solve(const double *l, double *d, double *u, double *f, int n);

/**
 * Prepare the partial LU factorisation of a cyclic, tridiagonal,
 * diagonally-dominant, square matrix A.
 *
 * We actually store the factorisation of the tridiagonal matrix B, formed by
 * removing the periodic elements from the corners and setting the first and
 * last elements of the main diagonal to 2*d[0] and d[n-1] + u[n-1]*l[0]/d[0]
 * respectively.
 *
 * We define the vector g as [-d[0], 0, 0, ..., 0, u[n-1]].
 *
 * @param l lower diagonal, overwritten with lower diagonal of L
 * @param d main diagonal, overwritten with main diagonal of L
 * @param u upper diagonal, overwritten with upper diagonal of U
 * @param q overwritten with B \ g
 * @param n size of the matrix
 */
void cyclic_tri_lu_factorise(
    const double *l, double *d, double *u, double *q, int n
);

/**
 * Given a partial LU factorisation of a cyclic, tridiagonal, square matrix A
 * = LU, solves Ax = f in place.
 *
 * See `cyclic_tri_lu_factorise` for the format of the matrix factorisation.
 *
 * @param l lower diagonal of L
 * @param d main diagonal of L
 * @param u upper diagonal of U
 * @param q B \ g
 * @param f right-hand side vector, overwritten with the solution
 * @param n size of the matrix
 */
void cyclic_tri_lu_solve(
    const double *l, const double *d, const double *u, const double *q,
    double *f, int n
);

/**
 * Given a cyclic, tridiagonal, square matrix A = LU, solves Ax = f in place.
 * Stores the partial LU factorisation of A in place so that it can be reused.
 *
 * See `cyclic_tri_lu_factorise` for the format of the matrix factorisation.
 *
 * @param l lower diagonal, overwritten with lower diagonal of L
 * @param d main diagonal, overwritten with main diagonal of L
 * @param u upper diagonal, overwritten with upper diagonal of U
 * @param q overwritten with B \ g
 * @param f right-hand side vector, overwritten with the solution
 * @param n size of the matrix
 */
void cyclic_tri_solve(
    const double *l, double *d, double *u, double *q, double *f, int n
);

#endif // TRI_SOLVE_H
