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

  END_TEST();
}
