// Included from tiny_plain.h
// DO NOT INCLUDE THIS FILE DIRECTLY.

#define SCITBX_ARRAY_FAMILY_TINY_CONVENIENCE_CONSTRUCTORS(class_name) \
explicit \
class_name( \
  value_type const& v0 \
) { \
  this->elems[0] = v0; \
} \
class_name( \
  value_type const& v0, \
  value_type const& v1 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
} \
class_name( \
  value_type const& v0, \
  value_type const& v1, \
  value_type const& v2 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
} \
class_name( \
  value_type const& v0, \
  value_type const& v1, \
  value_type const& v2, \
  value_type const& v3 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
} \
class_name( \
  value_type const& v0, \
  value_type const& v1, \
  value_type const& v2, \
  value_type const& v3, \
  value_type const& v4 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
} \
class_name( \
  value_type const& v0, \
  value_type const& v1, \
  value_type const& v2, \
  value_type const& v3, \
  value_type const& v4, \
  value_type const& v5 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
  this->elems[5] = v5; \
} \
class_name( \
  value_type const& v0, \
  value_type const& v1, \
  value_type const& v2, \
  value_type const& v3, \
  value_type const& v4, \
  value_type const& v5, \
  value_type const& v6 \
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
  value_type const& v0, \
  value_type const& v1, \
  value_type const& v2, \
  value_type const& v3, \
  value_type const& v4, \
  value_type const& v5, \
  value_type const& v6, \
  value_type const& v7 \
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
  value_type const& v0, \
  value_type const& v1, \
  value_type const& v2, \
  value_type const& v3, \
  value_type const& v4, \
  value_type const& v5, \
  value_type const& v6, \
  value_type const& v7, \
  value_type const& v8 \
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
  value_type const& v0, \
  value_type const& v1, \
  value_type const& v2, \
  value_type const& v3, \
  value_type const& v4, \
  value_type const& v5, \
  value_type const& v6, \
  value_type const& v7, \
  value_type const& v8, \
  value_type const& v9 \
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

#define SCITBX_ARRAY_FAMILY_TINY_COPY_AND_ASSIGNMENT(class_name) \
template <typename OtherElementType> \
class_name(tiny_plain<OtherElementType,N> const& rhs) { \
  copy_typeconv(rhs.begin(), rhs.end(), this->begin()); \
} \
template <typename OtherElementType> \
class_name<ElementType,N>& \
operator=(tiny_plain<OtherElementType,N> const& rhs) { \
  copy_typeconv(rhs.begin(), rhs.end(), this->begin()); \
  return *this; \
}
