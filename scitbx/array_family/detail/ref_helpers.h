// Included from ref.h
// DO NOT INCLUDE THIS FILE DIRECTLY.

#include <boost/operators.hpp>

namespace scitbx { namespace af { namespace detail {

  template<class T>
  class ref_reverse_iterator
    : public boost::random_access_iterator_helper<ref_reverse_iterator<T>, T>
  {
      T* p;
      typedef ref_reverse_iterator self;
      typedef std::ptrdiff_t Distance;
    public:
      ref_reverse_iterator() : p(0) {}
      explicit ref_reverse_iterator(T* beg) : p(beg-1) {}
      ref_reverse_iterator(self const& i) : p(i.p) {}
      self& operator=(self const& x) { p = x.p; return *this; }
      T& operator*() const { return *p; }
      bool operator==(self const& i) const { return p == i.p; }
      bool operator<(self const& i) const { return p > i.p; }
      self& operator++() { p--; return *this; }
      self& operator--() { p++; return *this; }
      self& operator+=(Distance n) { p -= n; return *this; }
      self& operator-=(Distance n) { p += n; return *this; }
      friend Distance operator-(self const& i, self const& j) {
        return j.p - i.p;
      }
  };
}}}

#define SCITBX_ARRAY_FAMILY_TYPEDEFS \
typedef ElementType                                     value_type; \
typedef ElementType*                                    iterator; \
typedef const ElementType*                              const_iterator; \
typedef detail::ref_reverse_iterator<ElementType>       reverse_iterator; \
typedef detail::ref_reverse_iterator<const ElementType> const_reverse_iterator; \
typedef ElementType&                                    reference; \
typedef ElementType const&                              const_reference; \
typedef std::size_t                                     size_type; \
typedef std::ptrdiff_t                                  difference_type;

#define SCITBX_ARRAY_FAMILY_BEGIN_END_ETC(this_type, beg, sz) \
ElementType* begin() { return beg; } \
const ElementType* begin() const { return beg; } \
ElementType* end() { return beg+sz; } \
const ElementType* end() const { return beg+sz; } \
reverse_iterator rbegin() { return reverse_iterator(beg+sz); } \
const_reverse_iterator rbegin() const { return const_reverse_iterator(beg+sz); } \
reverse_iterator rend() { return reverse_iterator(beg); } \
const_reverse_iterator rend() const { return const_reverse_iterator(beg); } \
ElementType& front() { return beg[0]; } \
ElementType const& front() const { return beg[0]; } \
ElementType& back() { return beg[sz-1]; } \
ElementType const& back() const { return beg[sz-1]; } \
 \
ElementType& operator[](size_type i) { return beg[i]; } \
ElementType const& operator[](size_type i) const { return beg[i]; } \
 \
ElementType& at(size_type i) { \
  if (i >= sz) throw_range_error(); return beg[i]; \
} \
ElementType const& at(size_type i) const { \
  if (i >= sz) throw_range_error(); return beg[i]; \
} \
this_type& fill(ElementType const& x) { \
  std::fill(this->begin(), this->end(), x); \
  return *this; \
}

#define SCITBX_ARRAY_FAMILY_TAKE_REF(beg, sz) \
af::ref<ElementType> \
ref() { \
  return af::ref<ElementType>(beg, sz); \
} \
af::const_ref<ElementType> \
const_ref() const { \
  return af::const_ref<ElementType>(beg, sz); \
}

#define SCITBX_ARRAY_FAMILY_TAKE_VERSA_REF(beg, ac) \
af::ref<ElementType, AccessorType> \
ref() { \
  return af::ref<ElementType, AccessorType>(beg, ac); \
} \
af::const_ref<ElementType, AccessorType> \
const_ref() const { \
  return af::const_ref<ElementType, AccessorType>(beg, ac); \
}
