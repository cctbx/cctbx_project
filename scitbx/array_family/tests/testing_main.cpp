// Simple main() function for cases where googletest is missing

#include "testing.h"

#ifdef NO_GTEST

int main(int argc, char const *argv[]) {
  return mockgtest::GetUnitTestImpl()->Run();
}

#endif