#ifndef TESTING_H
#define TESTING_H

#include <stdio.h>
#include <math.h>


const char *TESTN;
const char *SUBTESTN;
int ERR_COUNT;

#define START_TEST(title) { TESTN = title; SUBTESTN = ""; ERR_COUNT = 0; }


#define SUBTEST(subtitle) { SUBTESTN = subtitle; }


#define REQUIRE(cond) { int l = __LINE__;                                     \
    if (!(cond)) {                                                            \
      ERR_COUNT++;                                                            \
      fprintf(stderr, "  %s:%s (line %d) FAILED\n", TESTN, SUBTESTN, l);      \
    }                                                                         \
  }

// as above but exits early if the condition is not met
#define REQUIRE_BARRIER(cond) { int l = __LINE__;                             \
    if (!(cond)) {                                                            \
      ERR_COUNT++;                                                            \
      fprintf(stderr, "  %s:%s (line %d) FAILED\n", TESTN, SUBTESTN, l);      \
      fprintf(stderr, "    ABORTED EARLY DUE TO FAILURE\n");                  \
      return 1;                                                               \
    }                                                                         \
  }

#define REQUIRE_CLOSE(x, y, tol) { int l = __LINE__;                          \
    if (fabs(((double)(x)) - ((double)(y))) > tol) {                          \
      ERR_COUNT++;                                                            \
      fprintf(stderr, "  %s:%s (line %d) FAILED\n", TESTN, SUBTESTN, l);      \
    }                                                                         \
  }


#define END_TEST() { return ERR_COUNT; }


#endif //TESTING_H
