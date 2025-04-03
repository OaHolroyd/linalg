#include "testing.h"

#include "src/alloc.h"

int main(void) {
  START_TEST("alloc");

  /* check 2D allocation does not leak */
  SUBTEST("correct alloc") {
    double **A = malloc_d2d(10, 100);
    free_2d(A);
    // this should not cause a memory leak
  }

  /* check 2D allocation does not leak */
  SUBTEST("2D indexing") {
    // set up a 2D array with some values in
    int n = 10;
    int m = 3;
    double **A = malloc_d2d(n, m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        A[i][j] = i * m + j;
      }
    }

    // ensure that the values are correct regardless of the indexing
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        int k = i * m + j;
        REQUIRE_CLOSE(A[0][k], A[i][j], 1e-10);
      }
    }

    free_2d(A);
  }

  /* check calloc initialises to zero */
  SUBTEST("calloc") {
    int n = 10;
    int m = 3;
    double **A = calloc_d2d(n, m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        REQUIRE_BARRIER(A[i][j] == 0.0);
      }
    }
    free_2d(A);
  }

  END_TEST();
}
