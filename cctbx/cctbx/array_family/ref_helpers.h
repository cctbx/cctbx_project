// $Id$
// Included from ref.h
// DO NOT INCLUDE THIS FILE DIRECTLY.

#define CCTBX_ARRAY_FAMILY_TYPEDEFS \
typedef ElementType        value_type; \
typedef ElementType*       iterator; \
typedef const ElementType* const_iterator; \
typedef ElementType&       reference; \
typedef const ElementType& const_reference; \
typedef std::size_t        size_type; \
typedef std::ptrdiff_t     difference_type;

#define CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(beg, sz) \
ElementType* begin() { return beg; } \
const ElementType* begin() const { return beg; } \
ElementType* end() { return beg+sz; } \
const ElementType* end() const { return beg+sz; } \
ElementType& front() { return beg[0]; } \
const ElementType& front() const { return beg[0]; } \
ElementType& back() { return beg[sz-1]; } \
const ElementType& back() const { return beg[sz-1]; } \
 \
ElementType& operator[](size_type i) { return beg[i]; } \
const ElementType& operator[](size_type i) const { return beg[i]; } \
 \
ElementType& at(size_type i) { \
  if (i >= sz) throw_range_error(); return beg[i]; \
} \
const ElementType& at(size_type i) const { \
  if (i >= sz) throw_range_error(); return beg[i]; \
}

#define CCTBX_ARRAY_FAMILY_TAKE_REF(beg, sz) \
ref<ElementType> \
take_ref() { \
  return ref<ElementType>(this->elems, N); \
} \
const_ref<ElementType> \
take_const_ref() const { \
  return const_ref<ElementType>(this->elems, N); \
} \
operator ref<ElementType>() { \
  return ref<ElementType>(this->elems, N); \
} \
operator const_ref<ElementType>() const { \
  return const_ref<ElementType>(this->elems, N); \
}
