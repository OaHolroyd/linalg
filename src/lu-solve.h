#ifndef LU_SOLVE_H
#define LU_SOLVE_H

/**
 * Computes the LU factorisation of a matrix A with partial pivoting.
 *
 * This means that PA = LU, where P is a permutation matrix, L is unit lower
 * triangular and U is upper triangular. The matrix A is overwritten with the LU
 * factorisation (the overlapping of the diagonals of L and U is not a problem
 * since L has ones on the diagonal). This takes O(n^3) steps.
 *
 * @param A flattened matrix, overwritten with LU factorisation
 * @param piv pivot array, overwritten with pivot indices
 * @param n size of the matrix
 * @return 0 on success, row+1 on factorisation failure, -1 on other error
 */
int lu_factorise(double *A, int *piv, int n);

/**
 * Solves the system of equations LUx = Pf.
 *
 * Once the LU factorisation has been computed, the system can be solved in
 * O(n^2) steps:
 *  1. pivot the right-hand side vector f -> Pf
 *  2. solve Ly = Pf by forward substitution
 *  3. solve Ux = y by back substitution
 *
 * @param LU flattened matrix, containing LU factorisation
 * @param piv pivot array, containing pivot indices
 * @param f right-hand side vector, overwritten with solution
 * @param n size of the matrix
 * @return 0 on success, -1 on error
 */
int lu_solve_factorised(const double *LU, const int *piv, double *f, int n);

/**
 * Solves the system of equations Ax = f using LU factorisation with partial
 * pivoting.
 *
 * @param A flattened matrix, overwritten with LU factorisation
 * @param f right-hand side vector, overwritten with solution
 * @param piv pivot array, overwritten with pivot indices
 * @param n size of the matrix
 * @return 0 on success, -1 on error
 */
int lu_solve(double *A, double *f, int *piv, int n);

#endif // LU_SOLVE_H
