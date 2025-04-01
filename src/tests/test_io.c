#include "testing.h"

#include <math.h>
#include <unistd.h>

#include "src/alloc.h"
#include "src/io.h"

int main(void) {
  START_TEST("io");

  const int n = 3;
  const int m = 2;
  double **A = malloc_d2d(n, m);
  double **B = malloc_d2d(n, m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      A[i][j] = i * m + j + 0.1;
    }
  }

  /* valid output */
  SUBTEST("valid output") {
    const char filename[] = "tests/test_output.txt";
    int err = mat_outputf(filename, "%5.1lf", A[0], n, m);
    REQUIRE(err == 0);
    unlink(filename);
  }

  /* invalid output */
  SUBTEST("invalid output") {
    const char filename[] = "tests/nonexistent_dir/test_output.txt";
    int err = mat_outputf(filename, "%5.1lf", A[0], n, m);
    REQUIRE(err == 1);
  }

  /* read output back */
  SUBTEST("output/intput") {
    // output to a file
    const char filename[] = "tests/test_output.txt";
    int err = mat_outputf(filename, "%5.1lf", A[0], n, m);
    REQUIRE_BARRIER(err == 0);

    // read it back
    err = mat_input(filename, B[0], n, m);
    REQUIRE_BARRIER(err == 0);

    // check that the values are the same
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        REQUIRE(A[i][j] == B[i][j]);
      }
    }

    unlink(filename);
  }

  /* check output truncation */
  SUBTEST("truncated output") {
    // output to a file with no decimal places
    const char filename[] = "tests/test_output.txt";
    int err = mat_outputf(filename, "%5.0lf", A[0], n, m);
    REQUIRE_BARRIER(err == 0);

    // read it back
    err = mat_input(filename, B[0], n, m);
    REQUIRE_BARRIER(err == 0);

    // check that the values are rounded down
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        REQUIRE(floor(A[i][j]) == B[i][j]);
      }
    }

    unlink(filename);
  }

  /* check valid input */
  SUBTEST("valid input") {
    // output data to a file in the correct format
    const char filename[] = "tests/test_output.txt";
    FILE *fp = fopen(filename, "w");
    REQUIRE_BARRIER(fp != NULL);
    fprintf(fp, "1 2\n3 4\n5 6\n");
    fclose(fp);

    // read it back
    int err = mat_input(filename, B[0], n, m);
    REQUIRE_BARRIER(err == 0);

    // check that the values are as expected
    REQUIRE(B[0][0] == 1.0);
    REQUIRE(B[0][1] == 2.0);
    REQUIRE(B[1][0] == 3.0);
    REQUIRE(B[1][1] == 4.0);
    REQUIRE(B[2][0] == 5.0);
    REQUIRE(B[2][1] == 6.0);

    unlink(filename);
  }

  /* check valid input */
  SUBTEST("valid input (no newlines)") {
    // output data to a file in the correct format (even without newlines)
    const char filename[] = "tests/test_output.txt";
    FILE *fp = fopen(filename, "w");
    REQUIRE_BARRIER(fp != NULL);
    fprintf(fp, "1 2 3 4\n5 6\n");
    fclose(fp);

    // read it back
    int err = mat_input(filename, B[0], n, m);
    REQUIRE_BARRIER(err == 0);

    // check that the values are as expected
    REQUIRE(B[0][0] == 1.0);
    REQUIRE(B[0][1] == 2.0);
    REQUIRE(B[1][0] == 3.0);
    REQUIRE(B[1][1] == 4.0);
    REQUIRE(B[2][0] == 5.0);
    REQUIRE(B[2][1] == 6.0);

    unlink(filename);
  }

  /* check invalid inputs */
  SUBTEST("invalid input (too few entries)") {
    // output data to a file in an incorrect format (too few entries)
    const char filename[] = "tests/test_output.txt";
    FILE *fp = fopen(filename, "w");
    REQUIRE_BARRIER(fp != NULL);
    fprintf(fp, "1 2 3 4\n5\n");
    fclose(fp);

    // read it back
    int err = mat_input(filename, B[0], n, m);
    REQUIRE_BARRIER(err == 1);

    unlink(filename);
  }

  /* check invalid inputs */
  SUBTEST("invalid input (non-numeric data)") {
    // output data to a file in an incorrect format (too few entries)
    const char filename[] = "tests/test_output.txt";
    FILE *fp = fopen(filename, "w");
    REQUIRE_BARRIER(fp != NULL);
    fprintf(fp, "1 2 q3 4 5 6\n");
    fclose(fp);

    // read it back
    int err = mat_input(filename, B[0], n, m);
    REQUIRE_BARRIER(err == 1);

    unlink(filename);
  }

  free_2d(A);
  free_2d(B);

  END_TEST();
}
