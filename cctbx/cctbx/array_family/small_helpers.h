// $Id$
// Included from small_plain.h
// DO NOT INCLUDE THIS FILE DIRECTLY.

#define CCTBX_ARRAY_FAMILY_SMALL_CONSTRUCTORS(class_name) \
class_name() { this->m_size = 0; } \
explicit \
class_name(const size_type& sz, const ElementType& x = ElementType()) { \
  this->m_size = sz; \
  std::fill(this->begin(), this->end(), x); \
} \
template <typename OtherElementType> \
class_name(const OtherElementType* first, const OtherElementType* last) { \
  this->m_size = last - first; \
  copy_typeconv(first, last, this->begin()); \
}

#define CCTBX_ARRAY_FAMILY_SMALL_COPY_AND_ASSIGNMENT(class_name) \
template <typename OtherElementType, std::size_t OtherN> \
class_name(const small_plain<OtherElementType,OtherN>& rhs) { \
  this->m_size = rhs.size(); \
  if (this->m_size > N) throw_range_error(); \
  copy_typeconv(rhs.begin(), rhs.end(), this->begin()); \
} \
template <typename OtherElementType, std::size_t OtherN> \
class_name<ElementType,N>& \
operator=(const small_plain<OtherElementType,OtherN>& rhs) { \
  this->resize(rhs.size()); \
  copy_typeconv(rhs.begin(), rhs.end(), this->begin()); \
  return *this; \
}
