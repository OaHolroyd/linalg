#ifndef TESTING_H
#define TESTING_H

#include <math.h>
#include <stdio.h>

const char *TESTN;
const char *SUBTESTN;
int ERR_COUNT;

#define START_TEST(title)                                                      \
  {                                                                            \
    TESTN = title;                                                             \
    SUBTESTN = "";                                                             \
    ERR_COUNT = 0;                                                             \
  }

#define SUBTEST(subtitle)                                                      \
  {                                                                            \
    SUBTESTN = subtitle;                                                       \
    fprintf(stderr, "  %s:%s\n", TESTN, SUBTESTN);                             \
  }

#define REQUIRE(cond)                                                          \
  {                                                                            \
    if (!(cond)) {                                                             \
      ERR_COUNT++;                                                             \
      fprintf(stderr, "    line %d: FAILED\n", __LINE__);                      \
    }                                                                          \
  }

// as above but exits early if the condition is not met
#define REQUIRE_BARRIER(cond)                                                  \
  {                                                                            \
    if (!(cond)) {                                                             \
      ERR_COUNT++;                                                             \
      fprintf(stderr, "    line %d: ABORTED\n", __LINE__);                     \
      return 1;                                                                \
    }                                                                          \
  }

#define REQUIRE_CLOSE(x, y, tol)                                               \
  {                                                                            \
    if (fabs(((double)(x)) - ((double)(y))) > tol) {                           \
      ERR_COUNT++;                                                             \
      fprintf(stderr, "    line %d: FAILED\n", __LINE__);                      \
    }                                                                          \
  }

#define END_TEST()                                                             \
  {                                                                            \
    return ERR_COUNT;                                                          \
  }

#endif // TESTING_H
