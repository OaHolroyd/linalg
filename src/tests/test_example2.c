#include "testing.h"

#include "src/example.h"


int main(void) {
  START_TEST("example");

  /* check insertion works as expected */
  SUBTEST("example subtest") {
    REQUIRE(example_fn() == 1);
  }

  END_TEST();
}
