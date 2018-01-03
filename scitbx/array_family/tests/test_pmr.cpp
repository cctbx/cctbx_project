// Test the basic implementation of the C++17-style memory resources

#include <stdint.h>
#include <new>

#include <scitbx/array_family/polymorphic_allocator.h>

#include "testing.h"

using namespace scitbx::pmr;

TEST(PMRTest, BasicStatics) {
  // Test that new_delete_resource is unique
  ASSERT_NE(static_cast<void*>(NULL), new_delete_resource());
  ASSERT_EQ(new_delete_resource(), new_delete_resource());

  // Test that null_memory_resource is unique
  ASSERT_NE(static_cast<void*>(NULL), null_memory_resource());
  ASSERT_EQ(null_memory_resource(), null_memory_resource());

  // Ensure these aren't the same
  ASSERT_NE(null_memory_resource(), new_delete_resource());

  // Check that new_delete is the default
  ASSERT_EQ(get_default_resource(), new_delete_resource());
  // Try changing the default
  ASSERT_EQ(set_default_resource(null_memory_resource()),
            new_delete_resource());
  ASSERT_EQ(get_default_resource(), null_memory_resource());
  ASSERT_EQ(set_default_resource(new_delete_resource()),
            null_memory_resource());
}

TEST(PMRTest, new_delete_tests) {
  memory_resource& resource = *new_delete_resource();
  void* some_mem = resource.allocate(2048);
  // Make sure we got a value
  ASSERT_NE(static_cast<void*>(NULL), some_mem);
  // Validate the alignment of this. Hard-code the default system alignement
  // (equivalent to alignof(max_align_t) )
  ASSERT_EQ(reinterpret_cast<uintptr_t>(some_mem) % 2 * sizeof(std::size_t), 0);
  resource.deallocate(some_mem, 2048);
}

TEST(PMRTest, null_memory_tests) {
  memory_resource& resource = *null_memory_resource();
  ASSERT_THROW(resource.allocate(1), std::bad_alloc);
  ASSERT_NO_THROW(resource.deallocate(reinterpret_cast<void*>(NULL), 0));
}

#ifdef NDEBUG
TEST(PMRTest, DISABLED_test_no_alignment) {
#else
TEST(PMRTest, test_no_alignment) {
#endif
  // Check for non-available allocation if we're in debug mode
  memory_resource& resource = *new_delete_resource();
  ASSERT_DEATH(resource.allocate(100, 2), "Assertion failed");
}
