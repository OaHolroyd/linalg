#include "testing.h"

#include <unistd.h>

#include "src/alloc.h"
#include "src/io.h"

int main(void) {
  START_TEST("io");

  double **A = malloc_d2d(3, 2);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      A[i][j] = i * 2 + j + 0.1;
    }
  }

  /* valid output */
  SUBTEST("valid output") {
    const char filename[] = "tests/test_output.txt";
    int err = mat_outputf(filename, "%5.1lf", A[0], 3, 2);
    REQUIRE(err == 0);
    unlink(filename);
  }

  /* invalid output */
  SUBTEST("valid output") {
    const char filename[] = "tests/nonexistent_dir/test_output.txt";
    int err = mat_outputf(filename, "%5.1lf", A[0], 3, 2);
    REQUIRE(err == 1);
  }

  free_2d(A);

  END_TEST();
}
