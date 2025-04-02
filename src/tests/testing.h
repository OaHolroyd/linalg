#ifndef TESTING_H
#define TESTING_H

#include <stdio.h>
#include <math.h>


const char *TESTN;
const char *SUBTESTN;
int ERR_COUNT;

#define START_TEST(title) { TESTN = title; SUBTESTN = ""; ERR_COUNT = 0; }


#define SUBTEST(subtitle) {                      \
    SUBTESTN = subtitle;                       \
    fprintf(stderr, "  %s:%s\n", TESTN, SUBTESTN); \
  }


#define REQUIRE(cond) { int l = __LINE__;          \
    if (!(cond)) {                                 \
      ERR_COUNT++;                                 \
      fprintf(stderr, "    line %d: FAILED\n", l); \
    }                                              \
  }

// as above but exits early if the condition is not met
#define REQUIRE_BARRIER(cond) { int l = __LINE__;   \
    if (!(cond)) {                                  \
      ERR_COUNT++;                                  \
      fprintf(stderr, "    line %d: ABORTED\n", l); \
      return 1;                                     \
    }                                               \
  }

#define REQUIRE_CLOSE(x, y, tol) { int l = __LINE__; \
    if (fabs(((double)(x)) - ((double)(y))) > tol) { \
      ERR_COUNT++;                                   \
      fprintf(stderr, "    line %d: FAILED\n", l);   \
    }                                                \
  }


#define END_TEST() { return ERR_COUNT; }


#endif //TESTING_H
