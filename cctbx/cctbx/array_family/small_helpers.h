// $Id$
// Included from small_plain.h
// DO NOT INCLUDE THIS FILE DIRECTLY.

#define CCTBX_ARRAY_FAMILY_SMALL_CONSTRUCTORS(class_name) \
class_name() { this->m_size = 0; } \
explicit \
class_name(size_type n) { this->m_size = n; } \
class_name(size_type n, const ElementType& x) { \
  this->m_size = n; \
  std::fill(this->begin(), this->end(), x); \
}

#define CCTBX_ARRAY_FAMILY_SMALL_COPY_AND_ASSIGNMENT(class_name) \
template <typename OtherElementType, std::size_t OtherN> \
class_name(const small_plain<OtherElementType,OtherN>& rhs) { \
  this->m_size = rhs.size(); \
  if (this->m_size > N) throw_range_error(); \
  copy_typeconv(rhs.begin(), rhs.end(), this->elems); \
} \
template <typename OtherElementType, std::size_t OtherN> \
class_name<ElementType,N>& \
operator=(const small_plain<OtherElementType,OtherN>& rhs) { \
  this->resize(rhs.size()); \
  copy_typeconv(rhs.begin(), rhs.end(), this->elems); \
  return *this; \
}
