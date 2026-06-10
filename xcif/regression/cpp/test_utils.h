#ifndef XCIF_TEST_UTILS_H
#define XCIF_TEST_UTILS_H

#include <cstdio>

static int xcif_test_failures = 0;

#define CHECK(cond) \
  do { \
    if (!(cond)) { \
      std::fprintf(stderr, "FAIL: %s:%d: %s\n", __FILE__, __LINE__, #cond); \
      ++xcif_test_failures; \
    } \
  } while (0)

#define CHECK_EQ(a, b) CHECK((a) == (b))

#define CHECK_NE(a, b) CHECK((a) != (b))

#define CHECK_LT(a, b) CHECK((a) < (b))

#define CHECK_LE(a, b) CHECK((a) <= (b))

#define CHECK_GT(a, b) CHECK((a) > (b))

#define CHECK_GE(a, b) CHECK((a) >= (b))

#define XCIF_TEST_MAIN() \
  int main() { \
    run_all_tests(); \
    if (xcif_test_failures == 0) { \
      std::printf("All tests passed.\n"); \
    } else { \
      std::fprintf(stderr, "%d test(s) FAILED.\n", xcif_test_failures); \
    } \
    return xcif_test_failures > 0 ? 1 : 0; \
  }

#endif /* XCIF_TEST_UTILS_H */
