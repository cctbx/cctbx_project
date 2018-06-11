
#include <cstddef>
#include <scitbx/array_family/polymorphic_allocator.h>
#include <iostream>
#include <scitbx/array_family/shared_plain.h>
#include <vector>

#include "tst_af_helpers.cpp"

using namespace scitbx::pmr;
namespace af = scitbx::af;

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
  counting_resource(const counting_resource&);
  counting_resource& operator=(const counting_resource& rhs);

 public:
  counting_resource(memory_resource* parent = get_default_resource())
      : m_parent(parent), _allocated(0), _count(0), _deallocated(0) {}

  std::size_t count() { return _count; }
  std::size_t allocated() { return _allocated; }
  std::ptrdiff_t leaked() { return _allocated - _deallocated; }

  void reset() {
    _count = 0;
    _allocated = 0;
    _deallocated = 0;
  }
};

struct tst_test {
  static void exercise_1(void) {
    counting_resource resource;
    polymorphic_allocator<uint8_t> alloc(&resource);
    check_true(__LINE__, resource.count() ==  0);
    uint8_t* data = alloc.allocate(1024);
    check_true(__LINE__, resource.count() ==  1);
    check_true(__LINE__, resource.allocated() ==  1024);
    check_true(__LINE__, resource.leaked() ==  1024);
    alloc.deallocate(data, 1024);
    check_true(__LINE__, resource.leaked() ==  0);
    data = alloc.allocate(1024);
    resource.reset();
    check_true(__LINE__, resource.count() ==  0);
    check_true(__LINE__, resource.allocated() ==  0);
    check_true(__LINE__, resource.leaked() ==  0);
    alloc.deallocate(data, 1024);
    counting_resource child(&resource);
    resource.reset();
    void* vdata = child.allocate(1024);
    check_true(__LINE__, child.count() ==  resource.count());
    check_true(__LINE__, child.leaked() ==  resource.leaked());
    check_false(__LINE__, child.leaked() ==  0);
    child.deallocate(vdata, 1024);
    check_true(__LINE__, child.count() ==  resource.count());
  }

  static void exercise_2(void) {
    counting_resource resource;
    scitbx::af::sharing_handle a(&resource);
    check_true(__LINE__, resource.count() ==  0);
    check_true(__LINE__, *a.get_allocator().resource() ==  resource);
    scitbx::af::sharing_handle b(1e6, &resource);
    check_true(__LINE__, resource.count() ==  1);
    check_true(__LINE__, resource.allocated() ==  1e6);
    b.deallocate();
    check_false(__LINE__, resource.leaked());
    resource.reset();
    scitbx::af::sharing_handle* pb =
        new scitbx::af::sharing_handle(1e6, &resource);
    check_true(__LINE__, resource.allocated() ==  1e6);
    check_true(__LINE__, resource.leaked() ==  1e6);
    delete pb;
    check_false(__LINE__, resource.leaked());
    check_true(__LINE__, resource.count() ==  1);
  }

  static void exercise_6(void) {
    counting_resource resource;
    af::shared_plain<int>* a = new af::shared_plain<int>(&resource);
    check_true(__LINE__, resource.count() ==  1);
    a->push_back(4);
    check_true(__LINE__, resource.count() ==  3);
    check_true(__LINE__, a->capacity() ==  1);
    a->push_back(5);
    check_true(__LINE__, a->capacity() ==  2);
    a->push_back(6);
    check_true(__LINE__, resource.count() ==  7);
    check_true(__LINE__, a->size() ==  3);
    check_true(__LINE__, a->capacity() ==  4);
    a->reserve(100);
    check_true(__LINE__, resource.leaked() ==  sizeof(af::sharing_handle) + sizeof(int) * 100);
    delete a;
    check_true(__LINE__, resource.leaked() ==  0);
  }

  static void exercise_4(void) {
    counting_resource parent;
    counting_resource resource(&parent);
    af::shared_plain<int>* a = new af::shared_plain<int>(&resource);
    check_true(__LINE__, resource.count() ==  1);
    check_true(__LINE__, resource.allocated() ==  sizeof(af::sharing_handle));
    delete a;
    check_false(__LINE__, resource.leaked());
    resource.reset();
    a = new af::shared_plain<int>(5, &resource);
    check_true(__LINE__, a->size() ==  5);
    check_true(__LINE__, resource.count() ==  2);
    check_true(__LINE__, resource.allocated() ==  sizeof(af::sharing_handle) + sizeof(int) * 5);
    delete a;
    check_false(__LINE__, resource.leaked());
    resource.reset();
    a = new af::shared_plain<int>(af::reserve(5), &resource);
    check_true(__LINE__, a->size() ==  0);
    check_true(__LINE__, resource.count() ==  2);
    check_true(__LINE__, resource.allocated() ==  sizeof(af::sharing_handle) + sizeof(int) * 5);
    delete a;
    check_false(__LINE__, resource.leaked());
    resource.reset();
    a = new af::shared_plain<int>(7, 5, &resource);
    check_true(__LINE__, resource.count() ==  2);
    check_true(__LINE__, resource.allocated() ==  sizeof(af::sharing_handle) + sizeof(int) * 7);
    for (int i = 0; i < 7; ++i) {
      check_true(__LINE__, (*a)[i] == 5);
    }
    delete a;
    check_false(__LINE__, resource.leaked());
    resource.reset();
    a = new af::shared_plain<int>(10, af::init_functor_null<int>(), &resource);
    check_true(__LINE__, resource.count() ==  2);
    check_true(__LINE__, resource.allocated() ==
              sizeof(af::sharing_handle) + sizeof(int) * 10);
    delete a;
    check_false(__LINE__, resource.leaked());
    std::vector<int> sample;
    sample.push_back(4);
    sample.push_back(8);
    sample.push_back(15);
    sample.push_back(16);
    resource.reset();
    a = new af::shared_plain<int>(&(*sample.begin()), &(*sample.end()),
                                  &resource);
    check_true(__LINE__, resource.count() ==  2);
    check_true(__LINE__, resource.allocated() ==
              sizeof(af::sharing_handle) + sizeof(int) * sample.size());
    delete a;
    check_false(__LINE__, resource.leaked());
    std::vector<double> sampleD;
    sampleD.push_back(4.5);
    sampleD.push_back(8.5);
    sampleD.push_back(15.5);
    sampleD.push_back(16.5);
    resource.reset();
    a = new af::shared_plain<int>(&(*sampleD.begin()), &(*sampleD.end()),
                                  &resource);
    check_false(__LINE__,
        sizeof(double) ==
        sizeof(int));
    check_true(__LINE__, resource.count() ==  2);
    check_true(__LINE__, resource.allocated() ==
              sizeof(af::sharing_handle) + sizeof(int) * sample.size());
    delete a;
    check_false(__LINE__, resource.leaked());
    check_false(__LINE__, parent.leaked());
  }

  static void run(void) {
    exercise_1();
    exercise_2();
    exercise_6();
  }

  static void exercise_5(void) {
    counting_resource parent;
    counting_resource resource(&parent);
    af::shared_plain<int>* a = NULL;
    af::shared_plain<int>* root = new af::shared_plain<int>(10, &resource);
    check_false(__LINE__, parent.leaked() ==  0);
    resource.reset();
    a = new af::shared_plain<int>(*root);
    check_true(__LINE__, a->handle()->data ==  root->handle()->data);
    check_true(__LINE__, resource.count() ==  0);
    delete a;
    a = new af::shared_plain<int>(*root, af::weak_ref_flag());
    check_true(__LINE__, resource.count() ==  0);
    check_true(__LINE__, a->handle()->data ==  root->handle()->data);
    delete a;
    resource.reset();
    a = new af::shared_plain<int>(5, &resource);
    check_true(__LINE__, resource.count() ==  2);
    check_false(__LINE__, resource.leaked() ==  0);
    *a = *root;
    check_true(__LINE__, resource.count() ==  2);
    check_true(__LINE__, resource.leaked() ==  0);
    check_false(__LINE__, parent.leaked() ==  0);
    std::ptrdiff_t alloc_before = parent.leaked();
    delete a;
    check_true(__LINE__, parent.leaked() ==  alloc_before);
    counting_resource deep_copy_resource;
    {
      af::shared_plain<int> deepcopy = root->deep_copy(&deep_copy_resource);
      check_true(__LINE__, parent.leaked() ==  deep_copy_resource.leaked());
    }
    check_true(__LINE__, deep_copy_resource.leaked() ==  0);
    delete root;
    check_true(__LINE__, parent.leaked() ==  0);
  }

};

int main(int argc, char**)
{
  while(argc >= 0) {
    tst_test::run();
    tst_test::exercise_4();
    tst_test::exercise_5();
    if (argc < 2) break;
  }

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
