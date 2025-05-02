#ifndef PENT_SOLVE_H
#define PENT_SOLVE_H

/**
 * Factorises a pentadiagonal, diagonally dominant, square matrix A into A = LU.
 *
 * This is a special case of the LU factorisation, where the matrix A has five
 * diagonals (two lower, one main, two upper) and is diagonally dominant (i.e.
 * the absolute value of the diagonal element is greater than the sum of the
 * absolute values of the other elements in the same row). The resulting matrix
 * L is empty except for the main diagonal and the first two lower diagonals,
 * and the matrix U is empty except for the main diagonal and the first two
 * upper diagonals. Since the diagonal elements of U are all 1, it is not
 * necessary to store them so the factorisation can be performed in-place. This
 * takes O(n) steps, a significant improvement over the general O(n^3) LU
 * factorisation.
 *
 * @param l2 second lower diagonal, overwritten with second lower diagonal of L
 * @param l1 first lower diagonal, overwritten with first lower diagonal of L
 * @param d0 main diagonal, overwritten with main diagonal of L
 * @param u1 first upper diagonal, overwritten with first upper diagonal of U
 * @param u2 second upper diagonal, overwritten with second upper diagonal of U
 * @param n size of the matrix
 */
void pent_lu_factorise(
    const double *l2, double *l1, double *d0, double *u1, double *u2, int n
);

/**
 * Given an LU factorisation of a pentadiagonal, square, matrix A = LU, solves
 * Ax = f in place.
 *
 * Unlike in the general case (which is O(n^2)), this can be done in O(n) steps.
 *
 * @param l2 second lower diagonal of L
 * @param l1 first lower diagonal of L
 * @param l0 main diagonal of L
 * @param u1 first upper diagonal of U
 * @param u2 second upper diagonal of U
 * @param f right-hand side vector, overwritten with the solution
 * @param n size of the matrix
 */
void pent_lu_solve(
    const double *l2, const double *l1, const double *l0, const double *u1,
    const double *u2, double *f, int n
);

/**
 * Solves the system Ax = f in place, where A is a pentadiagonal.
 *
 * The index into the diagonal represents the row:
 *    d0_0  u1_0  u2_0     0     0     0     0    0 ...
 *    l1_1  d0_1  u1_1  u2_1     0     0     0    0 ...
 *    l2_2  l1_2  d0_2  u1_2  u2_2     0     0    0 ...
 *       0  l2_2  l1_2  d0_2  u1_2  u2_2     0    0 ...
 *       0     0  l2_2  l1_2  d0_2  u1_2  u2_2    0 ...
 * and so l1[0], l2[0], l2[1], u1[n-1], u2[n-2], u2[n-1] are not accessed.
 *
 * @param l2 second lower diagonal, overwritten with second lower diagonal of L
 * @param l1 first lower diagonal, overwritten with first lower diagonal of L
 * @param d0 main diagonal, overwritten with main diagonal of L
 * @param u1 first upper diagonal, overwritten with first upper diagonal of U
 * @param u2 second upper diagonal, overwritten with second upper diagonal of U
 * @param f right-hand side vector, overwritten with the solution
 * @param n size of the matrix
 */
void pent_solve(
    const double *l2, double *l1, double *d0, double *u1, double *u2, double *f,
    int n
);

/**
 * Prepare the partial LU factorisation of a cyclic, pentadiagonal,
 * diagonally-dominant, square matrix A.
 *
 * To perform the factorisation, the matrix is divided into blocks:
 *   A = [E  K]
 *       [H  C]
 * where
 *   E is a pentadiagonal n-2 by n-2 matrix,
 *   K is an n-2 by 2 matrix,
 *   H is an 2 by n-2 matrix,
 *   C is a 2 by 2 matrix.
 *
 * We LU factorise E and store E\K in K = [k0 k1]. In the leftover diagonal
 * elements no left unused by splitting off K, we store the 2 by 2 matrix C - H
 * E^-1 K.
 *
 * @param l2 second lower diagonal, overwritten with second lower diagonal of L
 * @param l1 first lower diagonal, overwritten with first lower diagonal of L
 * @param d0 main diagonal, overwritten with main diagonal of L
 * @param u1 first upper diagonal, overwritten with first upper diagonal of U
 * @param u2 second upper diagonal, overwritten with second upper diagonal of U
 * @param k0 overwritten with first column of E^-1 K
 * @param k1 overwritten with second column of E^-1 K
 * @param n size of the matrix
 */
void cyclic_pent_lu_factorise(
    const double *l2, double *l1, double *d0, double *u1, double *u2,
    double *k0, double *k1, int n
);

/**
 * Given a partial LU factorisation of a cyclic, pentadiagonal, square matrix A
 * = LU, solves Ax = f in place.
 *
 * See `cyclic_pent_lu_factorise` for the format of the matrix factorisation.
 *
 * @param l2 second lower diagonal of L
 * @param l1 first lower diagonal of L
 * @param l0 main diagonal of L
 * @param u1 first upper diagonal of U
 * @param u2 second upper diagonal of U
 * @param k0 the first column of E^-1 K
 * @param k1 the second column of E^-1 K
 * @param f right-hand side vector, overwritten with the solution
 * @param n size of the matrix
 */
void cyclic_pent_lu_solve(
    const double *l2, const double *l1, const double *l0, const double *u1,
    const double *u2, const double *k0, const double *k1, double *f, int n
);

/**
 * Given a cyclic, pentadiagonal, square matrix A = LU, solves Ax = f in place.
 * Stores the partial LU factorisation of A in place so that it can be reused.
 *
 * See `cyclic_pent_lu_factorise` for the format of the matrix factorisation.
 *
 * @param l2 second lower diagonal, overwritten with second lower diagonal of L
 * @param l1 first lower diagonal, overwritten with first lower diagonal of L
 * @param d0 main diagonal, overwritten with main diagonal of L
 * @param u1 first upper diagonal, overwritten with first upper diagonal of U
 * @param u2 second upper diagonal, overwritten with second upper diagonal of U
 * @param k0 overwritten with first column of E^-1 K
 * @param k1 overwritten with second column of E^-1 K
 * @param f right-hand side vector, overwritten with the solution
 * @param n size of the matrix
 */
void cyclic_pent_solve(
    const double *l2, double *l1, double *d0, double *u1, double *u2,
    double *k0, double *k1, double *f, int n
);

#endif // PENT_SOLVE_H
