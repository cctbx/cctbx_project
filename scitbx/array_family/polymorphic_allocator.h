/**
 * Backport of parts of C++17's allocator model.
 *
 * Proposed in https://isocpp.org/files/papers/N3916.pdf.
 * Accepted as
 * http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/n4480.html
 *
 * This is a much better, simpler model than in C++03,11,14,
 * But complete implementation requires C++11 support, or lots of TMP.
 * Thus, this is supposed to only approximate the behaviour so that when
 * C++11 is a minimal requirement we can swap out for a backport with
 * minimal difficulty.
 *
 * Definite differences:
 * - Alignment is currently ignored, as C++11 introduced a lot of syntax
 *   for managing this well. An aligned memory_resource allocator may
 *   be introduced in the future to manage with this, or just use aligned
 *   malloc. (Note that this means default alignment, not arbitrary)
 * - Full C++17 API is not implemented as it requires perfect forwarding
 *   and variadic templates - e.g. C++11. Possible with lots of TMP-based
 *   hoop jumping, and implemented in boost::container::pmr, but we don't
 *   want to introduce a dependency on an extra compiled module (and the
 *   complexity of boost's implementation is too much for our simple use
 *   case given that C++11 support should be coming... soon...)
 */

#ifndef SCITBX_MEMORY_RESOURCE_H
#define SCITBX_MEMORY_RESOURCE_H

#include <cassert>
#include <stdint.h> // For absolute-sized types (<cstdint> after C++11)
#include <new>      // For std::bad_alloc

#include <boost/config.hpp>  // compatibility macros
#include <boost/type_traits.hpp>    // Safety around StaticStore
#include <boost/static_assert.hpp>  //

namespace scitbx {
namespace pmr {

class memory_resource;
bool operator==(const memory_resource& a,
                const memory_resource& b) BOOST_NOEXCEPT;
bool operator!=(const memory_resource& a,
                const memory_resource& b) BOOST_NOEXCEPT;

template <class Tp>
class polymorphic_allocator;

template <class T1, class T2>
bool operator==(const polymorphic_allocator<T1>& a,
                const polymorphic_allocator<T2>& b) BOOST_NOEXCEPT;

template <class T1, class T2>
bool operator!=(const polymorphic_allocator<T1>& a,
                const polymorphic_allocator<T2>& b) BOOST_NOEXCEPT;

memory_resource* new_delete_resource() BOOST_NOEXCEPT;
memory_resource* null_memory_resource() BOOST_NOEXCEPT;
memory_resource* set_default_resource(memory_resource* r) BOOST_NOEXCEPT;
memory_resource* get_default_resource() BOOST_NOEXCEPT;

///////////////////////////////////////////////////////////////////////////////
// Class definitions

/// Abstract interface to represent a memory resource.
/// Based on C++17's std::pmr::memory_resource
class memory_resource {
  // This maximum align should be alignof(max_align_t) but 16 is the
  // requirement for 64-bit platforms apparently, and we don't have
  // access to this value in C++98. Apparently gcc's malloc uses
  // 2*sizeof(size_t) so that's what we use here.
  static BOOST_CONSTEXPR_OR_CONST size_t max_align = 2 * sizeof(std::size_t);

 public:
  /// Destroy this memory resource
  virtual ~memory_resource() {}

  void* allocate(std::size_t bytes, std::size_t alignment = max_align) {
    // Only allow default alignment
    assert(alignment == max_align);
    return do_allocate(bytes, alignment);
  }
  void deallocate(void* p, std::size_t bytes,
                  std::size_t alignment = max_align) {
    // Only allow default alignment
    assert(alignment == max_align);
    return do_deallocate(p, bytes, alignment);
  }

  bool is_equal(const memory_resource& other) const {
    return do_is_equal(other);
  }

 private:
  virtual void* do_allocate(std::size_t bytes, std::size_t alignment) = 0;
  virtual void do_deallocate(void* p, std::size_t bytes,
                             std::size_t alignment) = 0;

  virtual bool do_is_equal(const memory_resource& other) const = 0;
};

inline bool operator==(const memory_resource& a,
                       const memory_resource& b) BOOST_NOEXCEPT {
  return &a == &b || a.is_equal(b);
}

inline bool operator!=(const memory_resource& a,
                       const memory_resource& b) BOOST_NOEXCEPT {
  return !(a == b);
}

template <class Tp>
class polymorphic_allocator {
 private:
  memory_resource* memory_rsrc;

 public:
  typedef Tp value_type;
  typedef Tp* pointer;

  memory_resource* resource() const { return memory_rsrc; }

  // Constructors
  polymorphic_allocator() BOOST_NOEXCEPT : memory_rsrc(get_default_resource()) {
  }
  polymorphic_allocator(memory_resource* r) : memory_rsrc(r) {}
  polymorphic_allocator(const polymorphic_allocator& other)
      : memory_rsrc(other.resource()) {}

  template <class U>
  polymorphic_allocator(const polymorphic_allocator<U>& other) BOOST_NOEXCEPT
      : memory_rsrc(other.resource()) {}

  // Member functions
  Tp* allocate(size_t n) {
    // Divergence - don't have associated alignment pre-C++11
    return static_cast<Tp*>(
        memory_rsrc->allocate(n * sizeof(Tp) /*, alignof(Tp)*/));
  }
  void deallocate(Tp* p, size_t n) {
    // Divergence - don't have associated alignment pre-C++11
    memory_rsrc->deallocate(p, n * sizeof(Tp) /*, alignof(Tp)*/);
  }

  // Diversion from the standard
  //
  // Construct functions rely on variadic template arguments and perfect
  // forwarding, so we could skip their definitions.
  //
  // We still may need construct/destroy functions in order to wrap the
  // logic behind placement new/destuctor calling. This becomes somewhat
  // difficult because it's supposed to pass through the allocator if the
  // target type supports it (see: std::uses_allocator) and this is more
  // TMP than we want to put in for our minimal implementation. On the
  // other hand, we want to allow as much forward compatibility as possible
  // and removing excess allocator instances is less difficult than
  // manually rewriting manual destructs.
  //
  // So, making a suite of contruct functions, but not doing automatic
  // allocator detection and propogation

  template <class T>
  void construct(T* p) {
    new (p) T();
  }

  // Members taking a single parameter
  template <class T, class U>
  void construct(T* p, U& arg) {
    new (p) T(arg);
  }

  template <class T, class U>
  void construct(T* p, const U& arg) {
    new (p) T(arg);
  }

  // Members taking two parameters
  template <class T, class U1, class U2>
  void construct(T* p, U1& arg1, U2& arg2) {
    new (p) T(arg1, arg2);
  }
  template <class T, class U1, class U2>
  void construct(T* p, const U1& arg1, U2& arg2) {
    new (p) T(arg1, arg2);
  }
  template <class T, class U1, class U2>
  void construct(T* p, U1& arg1, const U2& arg2) {
    new (p) T(arg1, arg2);
  }
  template <class T, class U1, class U2>
  void construct(T* p, const U1& arg1, const U2& arg2) {
    new (p) T(arg1, arg2);
  }

  // Members taking three parameters
  template <class T, class U1, class U2, class U3>
  void construct(T* p, U1& arg1, U2& arg2, U3& arg3) {
    new (p) T(arg1, arg2, arg3);
  }
  template <class T, class U1, class U2, class U3>
  void construct(T* p, const U1& arg1, U2& arg2, U3& arg3) {
    new (p) T(arg1, arg2, arg3);
  }
  template <class T, class U1, class U2, class U3>
  void construct(T* p, U1& arg1, const U2& arg2, U3& arg3) {
    new (p) T(arg1, arg2, arg3);
  }
  template <class T, class U1, class U2, class U3>
  void construct(T* p, U1& arg1, U2& arg2, const U3& arg3) {
    new (p) T(arg1, arg2, arg3);
  }
  template <class T, class U1, class U2, class U3>
  void construct(T* p, const U1& arg1, const U2& arg2, U3& arg3) {
    new (p) T(arg1, arg2, arg3);
  }
  template <class T, class U1, class U2, class U3>
  void construct(T* p, const U1& arg1, U2& arg2, const U3& arg3) {
    new (p) T(arg1, arg2, arg3);
  }
  template <class T, class U1, class U2, class U3>
  void construct(T* p, U1& arg1, const U2& arg2, const U3& arg3) {
    new (p) T(arg1, arg2, arg3);
  }
  template <class T, class U1, class U2, class U3>
  void construct(T* p, const U1& arg1, const U2& arg2, const U3& arg3) {
    new (p) T(arg1, arg2, arg3);
  }

  template <class T>
  void destroy(T* p) {
    p->~T();
  }

  polymorphic_allocator select_on_container_copy_construction() const {
    return polymorphic_allocator();
  }

  BOOST_DELETED_FUNCTION(
      polymorphic_allocator& operator=(const polymorphic_allocator& rhs));
};

template <class T1, class T2>
bool operator==(const polymorphic_allocator<T1>& a,
                const polymorphic_allocator<T2>& b) BOOST_NOEXCEPT {
  return *a.resource() == *b.resource();
}

template <class T1, class T2>
bool operator!=(const polymorphic_allocator<T1>& a,
                const polymorphic_allocator<T2>& b) BOOST_NOEXCEPT {
  return !(a == b);
}

}  // namespace pmr
}  // namespace scitbx

///////////////////////////////////////////////////////////////////////////////
// Definitions that could be separated out

namespace scitbx {
namespace pmr {
namespace imp {

// Store implementation details for the default memory resources
class new_delete_memory_resource_impl : public memory_resource {
  virtual void* do_allocate(std::size_t bytes, std::size_t alignment) {
    return new uint8_t[bytes];
  }
  virtual void do_deallocate(void* p, std::size_t bytes,
                             std::size_t alignment) {
    delete[] static_cast<uint8_t**>(p);
  }

  virtual bool do_is_equal(const memory_resource& other) const BOOST_NOEXCEPT {
    return this == &other;
  }
};

class null_memory_resource_impl : public memory_resource {
  virtual void* do_allocate(std::size_t bytes, std::size_t alignment) {
    throw std::bad_alloc();
  }
  virtual void do_deallocate(void* p, std::size_t bytes,
                             std::size_t alignment) {}

  virtual bool do_is_equal(const memory_resource& other) const BOOST_NOEXCEPT {
    return this == &other;
  }
};

// Since we don't use SHARED libraries in this project, instead redeclaring
// all code (mostly inline) in every single MODULE library, there is no
// shared location to store the (mandated) global instance variables -
// the variables that would normally be declared in libc++ (or platform
// variation thereof).
//
// Static template class definitions seem to bypass the ODR rules,
// linking weakly so that the first physically loaded symbol location
// gets priority.
//
// The effect of this is that we can 'fake' a global variable by abusing
// the rules for static template definitions. This is a HACK HACK HACK
// and is only done as a last resort.
template <typename T>
struct StaticStore {
  // This must only ever be used for memory_resource
  BOOST_STATIC_ASSERT(( boost::is_same<T, memory_resource>::value ));

  static memory_resource* default_resource;
  static new_delete_memory_resource_impl* new_delete_resource;
  static null_memory_resource_impl* null_resource;
private:
  // This should never be instantiated
  StaticStore();
  StaticStore(const StaticStore&);
};

template <typename T>
new_delete_memory_resource_impl* StaticStore<T>::new_delete_resource =
    new new_delete_memory_resource_impl();
template <typename T>
memory_resource* StaticStore<T>::default_resource =
    StaticStore<T>::new_delete_resource;
template <typename T>
null_memory_resource_impl* StaticStore<T>::null_resource =
    new null_memory_resource_impl();

}  // namespace imp

inline memory_resource* new_delete_resource() BOOST_NOEXCEPT {
  return imp::StaticStore<memory_resource>::new_delete_resource;
}

inline memory_resource* null_memory_resource() BOOST_NOEXCEPT {
  return imp::StaticStore<memory_resource>::null_resource;
}

inline memory_resource* set_default_resource(memory_resource* r)
    BOOST_NOEXCEPT {
  memory_resource* old = imp::StaticStore<memory_resource>::default_resource;
  imp::StaticStore<memory_resource>::default_resource = r;
  return old;
}

inline memory_resource* get_default_resource() BOOST_NOEXCEPT {
  return imp::StaticStore<memory_resource>::default_resource;
}

}  // namespace pmr
}  // namespace scitbx

#endif
