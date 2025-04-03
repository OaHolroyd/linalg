#ifndef BLOCK_SOLVE_H
#define BLOCK_SOLVE_H

/**
 * Computes the solution to a system of linear equations decomposed into blocks.
 *
 * If R = [A B; C D], and f = [a; b], with A n x n, B n x m, C m x n, D m x m,
 * and a and b sized accordingly. Then we can use the method of Lu & Shiou 2000
 * to solve the system of equations Rx = f.
 *
 * This case assumes that both A and the Schur complement S = D - CA\B are
 * invertible.
 *
 * This solver uses a general LU factorisation of the blocks, and so will not be
 * any faster than just doing a full LU factorisation of the matrix R. However,
 * if A and/or S have some sort of structure (e.g. banded, triangular, etc.)
 * then the lu_solve calls can be replaced with more efficient versions, and the
 * block solve will be significantly faster.
 *
 * @param A upper left block
 * @param B upper right block
 * @param C lower left block
 * @param D lower right block
 * @param a upper right hand side
 * @param b lower right hand side
 * @param pivn pivot array for A
 * @param pivm pivot array for S
 * @param work interim workspace, should be at least size max(n*(m+1), (n+1)*m)
 * @param n upper left block size
 * @param m lower right block size
 * @return 0 on success, -1 on error
 */
int block_solve(
    double *A, const double *B, const double *C, double *D, double *a,
    double *b, int *pivn, int *pivm, double *work, int n, int m
);

/**
 * Example of how `block_solve` can be dramatically simplified for certain types
 * of blocks.
 *
 * In many cases (e.g. in solving PDEs), A = I and S can be computed directly
 *
 * @param B upper right block
 * @param C lower left block
 * @param S Schur complement, S = D - CA\B = D - CB
 * @param a upper right hand side
 * @param b lower right hand side
 * @param pivm pivot array for S
 * @param work interim workspace, should be at least size max(n, m)
 * @param n upper left block size
 * @param m lower right block size
 * @return 0 on success, -1 on error
 */
int block_solve_simplified(
    const double *B, const double *C, double *S, double *a, double *b,
    int *pivm, double *work, int n, int m
);

#endif // BLOCK_SOLVE_H
