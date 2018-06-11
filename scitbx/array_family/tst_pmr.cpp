#include <stdint.h>
#include <new>
#include <iostream>
#include <scitbx/array_family/polymorphic_allocator.h>
#include "tst_af_helpers.cpp"

using namespace scitbx::pmr;

struct array_exercise {
  static void run(void) {
    exercise_1();
  }
  static void exercise_1(void) {
    check_false(__LINE__, static_cast<void*>(NULL) == new_delete_resource());
    check_true(__LINE__, new_delete_resource() == new_delete_resource());
    check_false(__LINE__, static_cast<void*>(NULL) == null_memory_resource());
    check_true(__LINE__, null_memory_resource() == null_memory_resource());
    check_false(__LINE__, null_memory_resource() == new_delete_resource());
    check_true(__LINE__, get_default_resource() == new_delete_resource());
    check_true(__LINE__, set_default_resource(null_memory_resource()) == new_delete_resource());
    check_true(__LINE__, get_default_resource() == null_memory_resource());
    check_false(__LINE__, set_default_resource(new_delete_resource()) != null_memory_resource());
  }
};

struct array_exercise_3 {
  static void run() {
    exercise_2();
    exercise_2_2();
  }
  static void exercise_2(void) {
    memory_resource& resource = *new_delete_resource();
    void* some_mem = resource.allocate(2048);
    check_true(__LINE__, static_cast<void*>(NULL) != some_mem);
    check_true(__LINE__, reinterpret_cast<uintptr_t>(some_mem) % 2 * sizeof(std::size_t) == 0);
    resource.deallocate(some_mem, 2048);
  }

  static void exercise_2_2(void) {
    memory_resource& resource = *null_memory_resource();
    try {
      check_false(__LINE__, resource.allocate(1) || true);
    } catch (const std::bad_alloc &e) {
      check_true(__LINE__, true);
    }
    try {
      resource.deallocate(reinterpret_cast<void*>(NULL), 0);
      check_true(__LINE__, true);
    } catch (...) {
      check_false(__LINE__, true);
    }
  }
};

// #ifdef NDEBUG
// TEST(PMRTest, DISABLED_test_no_alignment) {
// #else
// TEST(PMRTest, test_no_alignment) {
// #endif
//   // Check for non-available allocation if we're in debug mode
//   memory_resource& resource = *new_delete_resource();
//   ASSERT_DEATH(resource.allocate(100, 2), "Assertion failed");
// }

int main(int argc, char* /*argv*/[])
{
  for(;;) {
    if (verbose) std::cout << __LINE__ << ":" << std::endl;
    array_exercise::run();
    if (verbose) std::cout << __LINE__ << ":" << std::endl;
    array_exercise_3::run();
    if (argc == 1) break;
  }

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
