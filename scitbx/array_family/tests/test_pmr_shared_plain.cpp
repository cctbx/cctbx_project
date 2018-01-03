// Test the addition of the polymorphic resources to array_family::shared_plain

#include <cstddef>  // size_t
#include <iostream>
#include <vector>

#include <boost/foreach.hpp>

#include <scitbx/array_family/polymorphic_allocator.h>
#include <scitbx/array_family/shared_plain.h>

#include "testing.h"

using namespace scitbx::pmr;
namespace af = scitbx::af;

// Simple test class to retain information about what was allocated
class counting_resource : public memory_resource {
 private:
  memory_resource* m_parent;
  std::size_t _allocated;
  std::size_t _count;
  std::size_t _deallocated;

  virtual void* do_allocate(std::size_t bytes, std::size_t alignment) {
    _count += 1;
    _allocated += bytes;
    return m_parent->allocate(bytes, alignment);
  }
  virtual void do_deallocate(void* p, std::size_t bytes,
                             std::size_t alignment) {
    _deallocated += bytes;
    m_parent->deallocate(p, bytes, alignment);
  }

  virtual bool do_is_equal(const memory_resource& other) const BOOST_NOEXCEPT {
    return this == &other;
  }

  // Don't allow copying of these
  counting_resource(const counting_resource&);
  counting_resource& operator=(const counting_resource& rhs);

 public:
  counting_resource(memory_resource* parent = get_default_resource())
      : m_parent(parent), _allocated(0), _count(0), _deallocated(0) {}

  /// Get the number of allocations made
  std::size_t count() { return _count; }
  /// Get the number of bytes allocated
  std::size_t allocated() { return _allocated; }
  /// Get the total memory currently unfreed
  std::ptrdiff_t leaked() { return _allocated - _deallocated; }

  void reset() {
    _count = 0;
    _allocated = 0;
    _deallocated = 0;
  }
};

// Test that our diagnostic tool works
TEST(PMR_shared_plain, TestReportingResource) {
  counting_resource resource;
  polymorphic_allocator<uint8_t> alloc(&resource);
  ASSERT_EQ(resource.count(), 0);
  uint8_t* data = alloc.allocate(1024);
  ASSERT_EQ(resource.count(), 1);
  ASSERT_EQ(resource.allocated(), 1024);
  ASSERT_EQ(resource.leaked(), 1024);
  alloc.deallocate(data, 1024);
  ASSERT_EQ(resource.leaked(), 0);

  // Check resetting
  data = alloc.allocate(1024);
  resource.reset();
  ASSERT_EQ(resource.count(), 0);
  ASSERT_EQ(resource.allocated(), 0);
  ASSERT_EQ(resource.leaked(), 0);
  alloc.deallocate(data, 1024);

  // Chaining
  counting_resource child(&resource);
  resource.reset();
  void* vdata = child.allocate(1024);
  ASSERT_EQ(child.count(), resource.count());
  ASSERT_EQ(child.leaked(), resource.leaked());
  ASSERT_NE(child.leaked(), 0);
  child.deallocate(vdata, 1024);
  ASSERT_EQ(child.count(), resource.count());
}

TEST(PMR_shared_plain, SharingHandleWithAllocator) {
  counting_resource resource;
  // Basic construction should not cause allocation
  scitbx::af::sharing_handle a(&resource);
  ASSERT_EQ(resource.count(), 0);
  // Basic consutrction SHOULD pass through the allocator
  ASSERT_EQ(*a.get_allocator().resource(), resource);

  // Give me ONE MILLION of your precious bytes
  scitbx::af::sharing_handle b(1e6, &resource);
  ASSERT_EQ(resource.count(), 1);
  ASSERT_EQ(resource.allocated(), 1e6);
  b.deallocate();
  ASSERT_FALSE(resource.leaked());

  // Test deconstruction deallocates
  resource.reset();
  scitbx::af::sharing_handle* pb =
      new scitbx::af::sharing_handle(1e6, &resource);
  ASSERT_EQ(resource.allocated(), 1e6);
  ASSERT_EQ(resource.leaked(), 1e6);
  delete pb;
  ASSERT_FALSE(resource.leaked());
  // Total allocations maded?
  ASSERT_EQ(resource.count(), 1);
}

TEST(PMR_shared_plain, constructors) {
  counting_resource parent;
  counting_resource resource(&parent);

  // Empty constructor
  af::shared_plain<int>* a = new af::shared_plain<int>(&resource);
  ASSERT_EQ(resource.count(), 1);
  ASSERT_EQ(resource.allocated(), sizeof(af::sharing_handle));
  delete a;
  ASSERT_FALSE(resource.leaked());

  // Size constructor
  resource.reset();
  a = new af::shared_plain<int>(5, &resource);
  ASSERT_EQ(a->size(), 5);
  ASSERT_EQ(resource.count(), 2);
  ASSERT_EQ(resource.allocated(), sizeof(af::sharing_handle) + sizeof(int) * 5);
  delete a;
  ASSERT_FALSE(resource.leaked());

  // Reserved size constructor
  resource.reset();
  a = new af::shared_plain<int>(af::reserve(5), &resource);
  ASSERT_EQ(a->size(), 0);
  ASSERT_EQ(resource.count(), 2);
  ASSERT_EQ(resource.allocated(), sizeof(af::sharing_handle) + sizeof(int) * 5);
  delete a;
  ASSERT_FALSE(resource.leaked());

  // Size and Initial Element Constructor
  resource.reset();
  a = new af::shared_plain<int>(7, 5, &resource);
  ASSERT_EQ(resource.count(), 2);
  ASSERT_EQ(resource.allocated(), sizeof(af::sharing_handle) + sizeof(int) * 7);
  // Make sure every element was set
  BOOST_FOREACH (int item, *a) { ASSERT_EQ(item, 5); }
  delete a;
  ASSERT_FALSE(resource.leaked());

  // - Size, Initialization functor Constructor
  resource.reset();
  a = new af::shared_plain<int>(10, af::init_functor_null<int>(), &resource);
  ASSERT_EQ(resource.count(), 2);
  ASSERT_EQ(resource.allocated(),
            sizeof(af::sharing_handle) + sizeof(int) * 10);
  delete a;
  ASSERT_FALSE(resource.leaked());

  // Iterator begin, end Constructor - note that uses pointers not iterator
  // references
  std::vector<int> sample;
  sample.push_back(4);
  sample.push_back(8);
  sample.push_back(15);
  sample.push_back(16);
  resource.reset();
  a = new af::shared_plain<int>(&(*sample.begin()), &(*sample.end()),
                                &resource);
  ASSERT_EQ(resource.count(), 2);
  ASSERT_EQ(resource.allocated(),
            sizeof(af::sharing_handle) + sizeof(int) * sample.size());
  delete a;
  ASSERT_FALSE(resource.leaked());

  // - Templated Iterator begin, end Conversion Constructor
  std::vector<double> sampleD;
  sampleD.push_back(4.5);
  sampleD.push_back(8.5);
  sampleD.push_back(15.5);
  sampleD.push_back(16.5);
  resource.reset();
  a = new af::shared_plain<int>(&(*sampleD.begin()), &(*sampleD.end()),
                                &resource);
  ASSERT_NE(
      sizeof(double),
      sizeof(int));  // This tests better validated if size change required
  ASSERT_EQ(resource.count(), 2);
  ASSERT_EQ(resource.allocated(),
            sizeof(af::sharing_handle) + sizeof(int) * sample.size());
  delete a;
  ASSERT_FALSE(resource.leaked());

  ASSERT_FALSE(parent.leaked());
}

TEST(PMR_shared_plain, weak_and_copy_constructors) {
  counting_resource parent;
  counting_resource resource(&parent);

  af::shared_plain<int>* a = NULL;

  // Create a base shared_plain to share references to
  af::shared_plain<int>* root = new af::shared_plain<int>(10, &resource);
  ASSERT_NE(parent.leaked(), 0);

  // NONSTD Copy constructor should copy nothing, no allocations made
  resource.reset();
  a = new af::shared_plain<int>(*root);
  ASSERT_EQ(a->handle()->data, root->handle()->data);
  ASSERT_EQ(resource.count(), 0);
  delete a;

  // NONSTD Weak copy constructor
  a = new af::shared_plain<int>(*root, af::weak_ref_flag());
  ASSERT_EQ(resource.count(), 0);
  ASSERT_EQ(a->handle()->data, root->handle()->data);
  delete a;

  // NONSTD Copy-assign (Share assign) operator
  resource.reset();
  a = new af::shared_plain<int>(5, &resource);
  ASSERT_EQ(resource.count(), 2);
  ASSERT_NE(resource.leaked(), 0);
  *a = *root;
  // Should have copied as a strong reference
  ASSERT_EQ(resource.count(), 2);
  // ...but returned everything we had allocated before
  ASSERT_EQ(resource.leaked(), 0);
  // Make sure that deleting this touches nothing
  ASSERT_NE(parent.leaked(), 0);
  std::ptrdiff_t alloc_before = parent.leaked();
  delete a;
  ASSERT_EQ(parent.leaked(), alloc_before);

  // Now, do a deep copy - equivalent of a normal copy
  // But do it to a new resource
  counting_resource deep_copy_resource;
  // Deep copy only returns an object directly so put it in a
  // sub-scope so that we can ensure it's storage duration closes
  {
    af::shared_plain<int> deepcopy = root->deep_copy(&deep_copy_resource);
    ASSERT_EQ(parent.leaked(), deep_copy_resource.leaked());
  }
  ASSERT_EQ(deep_copy_resource.leaked(), 0);
  delete root;
  ASSERT_EQ(parent.leaked(), 0);
}

TEST(PMR_shared_plain, resizing_methods_and_effects) {
  counting_resource resource;

  af::shared_plain<int>* a = new af::shared_plain<int>(&resource);
  ASSERT_EQ(resource.count(), 1);
  // Creates a new, sized array then swaps; so two allocs per increase
  a->push_back(4);
  ASSERT_EQ(resource.count(), 3);
  ASSERT_EQ(a->capacity(), 1);
  a->push_back(5);
  ASSERT_EQ(a->capacity(), 2);
  a->push_back(6);
  ASSERT_EQ(resource.count(), 7);
  ASSERT_EQ(a->size(), 3);
  ASSERT_EQ(a->capacity(), 4);

  a->reserve(100);
  ASSERT_EQ(resource.leaked(), sizeof(af::sharing_handle) + sizeof(int) * 100);
  delete a;
  ASSERT_EQ(resource.leaked(), 0);
}

// - array_adaptor<OtherArrayType> Copy Constructor
// - DEEP Copy
// - weak_ref copying