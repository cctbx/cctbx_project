// $Id$
// Included from tiny_plain.h
// DO NOT INCLUDE THIS FILE DIRECTLY.

#define CCTBX_ARRAY_FAMILY_TINY_CONVENIENCE_CONSTRUCTORS(class_name) \
explicit \
class_name( \
  const value_type& v0 \
) { \
  this->elems[0] = v0; \
} \
class_name( \
  const value_type& v0, \
  const value_type& v1 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
} \
class_name( \
  const value_type& v0, \
  const value_type& v1, \
  const value_type& v2 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
} \
class_name( \
  const value_type& v0, \
  const value_type& v1, \
  const value_type& v2, \
  const value_type& v3 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
} \
class_name( \
  const value_type& v0, \
  const value_type& v1, \
  const value_type& v2, \
  const value_type& v3, \
  const value_type& v4 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
} \
class_name( \
  const value_type& v0, \
  const value_type& v1, \
  const value_type& v2, \
  const value_type& v3, \
  const value_type& v4, \
  const value_type& v5 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
  this->elems[5] = v5; \
} \
class_name( \
  const value_type& v0, \
  const value_type& v1, \
  const value_type& v2, \
  const value_type& v3, \
  const value_type& v4, \
  const value_type& v5, \
  const value_type& v6 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
  this->elems[5] = v5; \
  this->elems[6] = v6; \
} \
class_name( \
  const value_type& v0, \
  const value_type& v1, \
  const value_type& v2, \
  const value_type& v3, \
  const value_type& v4, \
  const value_type& v5, \
  const value_type& v6, \
  const value_type& v7 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
  this->elems[5] = v5; \
  this->elems[6] = v6; \
  this->elems[7] = v7; \
} \
class_name( \
  const value_type& v0, \
  const value_type& v1, \
  const value_type& v2, \
  const value_type& v3, \
  const value_type& v4, \
  const value_type& v5, \
  const value_type& v6, \
  const value_type& v7, \
  const value_type& v8 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
  this->elems[5] = v5; \
  this->elems[6] = v6; \
  this->elems[7] = v7; \
  this->elems[8] = v8; \
} \
class_name( \
  const value_type& v0, \
  const value_type& v1, \
  const value_type& v2, \
  const value_type& v3, \
  const value_type& v4, \
  const value_type& v5, \
  const value_type& v6, \
  const value_type& v7, \
  const value_type& v8, \
  const value_type& v9 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
  this->elems[5] = v5; \
  this->elems[6] = v6; \
  this->elems[7] = v7; \
  this->elems[8] = v8; \
  this->elems[9] = v9; \
} \
template <typename OtherElementType> \
class_name(const OtherElementType* first, const OtherElementType* last) { \
  if (last - first != this->size()) throw_range_error(); \
  copy_typeconv(first, last, this->begin()); \
}

#define CCTBX_ARRAY_FAMILY_TINY_COPY_AND_ASSIGNMENT(class_name) \
template <typename OtherElementType> \
class_name(const tiny_plain<OtherElementType,N>& rhs) { \
  copy_typeconv(rhs.begin(), rhs.end(), this->begin()); \
} \
template <typename OtherElementType> \
class_name<ElementType,N>& \
operator=(const tiny_plain<OtherElementType,N>& rhs) { \
  copy_typeconv(rhs.begin(), rhs.end(), this->begin()); \
  return *this; \
}
