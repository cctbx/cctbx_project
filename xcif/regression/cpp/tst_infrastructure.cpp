/* Smoke test for the xcif C++ test infrastructure (test_utils.h macros).
   No xcif parser code is exercised here — this only verifies that the
   CHECK/CHECK_EQ/XCIF_TEST_MAIN macros compile, link, and report correctly. */
#include "test_utils.h"

static void test_check_passes() {
  CHECK(1 == 1);
  CHECK(true);
}

static void test_check_eq_passes() {
  int a = 42;
  int b = 42;
  CHECK_EQ(a, b);
}

static void test_comparisons() {
  CHECK_NE(1, 2);
  CHECK_LT(1, 2);
  CHECK_LE(1, 1);
  CHECK_GT(2, 1);
  CHECK_GE(2, 2);
}

static void run_all_tests() {
  test_check_passes();
  test_check_eq_passes();
  test_comparisons();
}

XCIF_TEST_MAIN()
