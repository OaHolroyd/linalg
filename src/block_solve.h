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
    double *A, double *B, const double *C, double *D, double *a, double *b,
    int *pivn, int *pivm, double *work, int n, int m
);

#endif // BLOCK_SOLVE_H
