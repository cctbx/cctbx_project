/**
 * Basic testing infrastructure for array_family polymorphic resources.
 *
 * Provides an extremely minimal replacement for some of the GoogleTest
 * testing API, under the hope that at some point a proper C++ testing
 * framework will be included with the distribution.
 *
 * If googletest is present, determined by the define NO_GTEST being
 # undefined, then use that instead of our minimal replacement.
*/

#ifndef TESTING_WRAPPER_H
#define TESTING_WRAPPER_H

// Unless we've been asked to not use GTest
#ifndef NO_GTEST
#include <gtest/gtest.h>
#else

#include <iostream>
#include <string>
#include <vector>

#include <boost/algorithm/string/predicate.hpp>

#include <scitbx/error.h>

namespace mockgtest {

// Base class for test bodies
class Test {
 public:
  virtual ~Test(){};
  virtual void TestBody(void) = 0;
};

// Simple polymorphic factory to make a test on demand
class TestFactoryBase {
 public:
  virtual Test* create() = 0;
};
template <class T>
class TestFactory : public TestFactoryBase {
 public:
  TestFactory(){};
  virtual Test* create() { return new T(); }
};

// Hold information about a single test
class TestInfo {
 public:
  TestInfo(std::string casename, std::string name, TestFactoryBase* factory)
      : testName(name), testCase(casename), factory(factory) {}
  std::string testName;
  std::string testCase;
  TestFactoryBase* factory;
};

// Store list of test instances
class UnitTestImpl {
  std::vector<TestInfo*> _tests;

 public:
  void AddTestInfo(TestInfo* test) { _tests.push_back(test); }

  int Run(void) {
    for (std::vector<TestInfo*>::iterator it = _tests.begin();
         it != _tests.end(); ++it) {
      if (boost::starts_with((*it)->testName, "DISABLED")) {
        std::cout << "(Skipping) " << (*it)->testCase << "::" << (*it)->testName
                  << std::endl;
      } else {
        std::cout << "Running    " << (*it)->testCase << "::" << (*it)->testName
                  << std::endl;
        Test* test = (*it)->factory->create();
        test->TestBody();
        delete test;
      }
    }
    return 0;
  }
};

// Hackily abuse template ODR to allow static stores without shared code
// libraries
template <typename T>
struct UnitTestStaticStore {
  static UnitTestImpl* unittestimpl;
};
template <typename T>
UnitTestImpl* UnitTestStaticStore<T>::unittestimpl = new UnitTestImpl();

// Get a static test implementation
inline UnitTestImpl* GetUnitTestImpl() {
  return UnitTestStaticStore<UnitTestImpl>::unittestimpl;
}

// Create and register a static unit case
inline TestInfo* MakeAndRegisterTestInfo(std::string casename, std::string name,
                                         TestFactoryBase* factory) {
  TestInfo* ti = new TestInfo(casename, name, factory);
  GetUnitTestImpl()->AddTestInfo(ti);
  return ti;
}

#define TEST_CLASS_NAME_(test_case_name, test_name) \
  test_case_name##_##test_name##_Test

#define TEST(test_case_name, test_name)                                        \
  class TEST_CLASS_NAME_(test_case_name, test_name) : public mockgtest::Test { \
   public:                                                                     \
    TEST_CLASS_NAME_(test_case_name, test_name)(){};                           \
    virtual void TestBody();                                                   \
    static mockgtest::TestInfo* const _test_info;                              \
  };                                                                           \
  mockgtest::TestInfo* const TEST_CLASS_NAME_(test_case_name,                  \
                                              test_name)::_test_info =         \
      mockgtest::MakeAndRegisterTestInfo(                                      \
          #test_case_name, #test_name,                                         \
          new mockgtest::TestFactory<TEST_CLASS_NAME_(test_case_name,          \
                                                      test_name)>);            \
  void TEST_CLASS_NAME_(test_case_name, test_name)::TestBody()

// void TEST_CLASS_NAME_(package, test_name)(void)

#define ASSERT_EQ(l, r) SCITBX_ASSERT(l == r)
#define ASSERT_NE(l, r) SCITBX_ASSERT(l != r)
#define ASSERT_FALSE(expr) SCITBX_ASSERT(!expr)

#define ASSERT_THROW(expr, expected_exception) \
  try {                                        \
    expr;                                      \
    SCITBX_ASSERT(0);                          \
  } catch (expected_exception & e) {           \
  } catch (...) {                              \
    SCITBX_ASSERT(0);                          \
  }

#define ASSERT_NO_THROW(expr) \
  try {                       \
    expr;                     \
  } catch (...) {             \
    SCITBX_ASSERT(0);         \
  }

}  // namespace mockgtest
#define ASSERT_DEATH(stmt, regex)                                      \
  std::cout << "Warning: Cannot run death test on " << __FILE__ << ":" \
            << __LINE__ << std::endl;

#endif

#endif  // inclusion guard