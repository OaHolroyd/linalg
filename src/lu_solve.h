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
 * @param piv pivot array, overwritten with pivot indices. If NULL, assumes no
 * pivoting.
 * @param n size of the matrix
 * @return 0 on success, row+1 on factorisation failure, -1 on other error
 */
int lu_factorise(double *A, int *piv, int n);

/**
 * Computes the LU factorisation of a matrix A with no pivoting.
 *
 * @param A flattened matrix, overwritten with LU factorisation
 * @param n size of the matrix
 * @return 0 on success, row+1 on factorisation failure, -1 on other error
 */
int lu_factorise_no_pivoting(double *A, int n);

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
 * @param piv pivot array, containing pivot indices. If NULL, assumes no
 * pivoting.
 * @param f right-hand side vector, overwritten with solution
 * @param n size of the matrix
 */
void lu_solve_factorised(const double *LU, const int *piv, double *f, int n);

/**
 * Solves the system of equations LUX = PF.
 *
 * As for `lu_solve_factorised` but for multiple right-hand side vectors.
 *
 * @param LU flattened matrix, containing LU factorisation
 * @param piv pivot array, containing pivot indices. If NULL, assumes no
 * pivoting.
 * @param F right-hand side vectors, overwritten with solution
 * @param n number of rows of the matrix
 * @param m number of right-hand side vectors
 */
void lu_solve_factorised_multi(
    const double *LU, const int *piv, double *F, int n, int m
);

/**
 * Solves the system of equations LUx = f.
 *
 * @param LU flattened matrix, containing LU factorisation
 * @param f right-hand side vector, overwritten with solution
 * @param n size of the matrix
 */
void lu_solve_factorised_no_pivoting(const double *LU, double *f, int n);

/**
 * Solves the system of equations Ax = f using LU factorisation with partial
 * pivoting.
 *
 * @param A flattened matrix, overwritten with LU factorisation
 * @param f right-hand side vector, overwritten with solution
 * @param piv pivot array, overwritten with pivot indices. If NULL, assumes no
 * pivoting.
 * @param n size of the matrix
 * @return 0 on success, -1 on error
 */
int lu_solve(double *A, double *f, int *piv, int n);

/**
 * Solves the system of equations AX = F using LU factorisation with partial
 * pivoting.
 *
 * @param A flattened matrix, overwritten with LU factorisation
 * @param F right-hand side vectors, overwritten with solution
 * @param piv pivot array, overwritten with pivot indices. If NULL, assumes no
 * pivoting.
 * @param n size of the matrix
 * @param m number of right-hand side vectors
 * @return 0 on success, -1 on error
 */
int lu_solve_multi(double *A, double *F, int *piv, int n, int m);

/**
 * Solves the system of equations Ax = f using LU factorisation with no
 * pivoting.
 *
 * @param A flattened matrix, overwritten with LU factorisation
 * @param f right-hand side vector, overwritten with solution
 * @param n size of the matrix
 * @return 0 on success, -1 on error
 */
int lu_solve_no_pivoting(double *A, double *f, int n);

#endif // LU_SOLVE_H
