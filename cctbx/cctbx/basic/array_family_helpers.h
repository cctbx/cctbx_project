// Included from array_family.h
// DO NOT INCLUDE THIS FILE DIRECTLY.

#define CCTBX_BASIC_ARRAY_FAMILY_TYPEDEFS \
typedef ElementType        value_type; \
typedef ElementType*       iterator; \
typedef const ElementType* const_iterator; \
typedef ElementType&       reference; \
typedef const ElementType& const_reference; \
typedef std::size_t        size_type; \
typedef std::ptrdiff_t     difference_type;

#define CCTBX_BASIC_ARRAY_FAMILY_CONVENIENCE_CONSTRUCTORS_FIXSIZE(class_name) \
explicit \
class_name( \
  const ElementType& v0 \
) { \
  this->elems[0] = v0; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3, \
  const ElementType& v4 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3, \
  const ElementType& v4, \
  const ElementType& v5 \
) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
  this->elems[5] = v5; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3, \
  const ElementType& v4, \
  const ElementType& v5, \
  const ElementType& v6 \
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
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3, \
  const ElementType& v4, \
  const ElementType& v5, \
  const ElementType& v6, \
  const ElementType& v7 \
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
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3, \
  const ElementType& v4, \
  const ElementType& v5, \
  const ElementType& v6, \
  const ElementType& v7, \
  const ElementType& v8 \
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
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3, \
  const ElementType& v4, \
  const ElementType& v5, \
  const ElementType& v6, \
  const ElementType& v7, \
  const ElementType& v8, \
  const ElementType& v9 \
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
}

#define CCTBX_BASIC_ARRAY_FAMILY_CONVENIENCE_CONSTRUCTORS_VARSIZE(class_name) \
explicit \
class_name( \
  const ElementType& v0 \
) : m_size(1) { \
  this->elems[0] = v0; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1 \
) : m_size(2) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2 \
) : m_size(3) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3 \
) : m_size(4) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3, \
  const ElementType& v4 \
) : m_size(5) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3, \
  const ElementType& v4, \
  const ElementType& v5 \
) : m_size(6) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
  this->elems[5] = v5; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3, \
  const ElementType& v4, \
  const ElementType& v5, \
  const ElementType& v6 \
) : m_size(7) { \
  this->elems[0] = v0; \
  this->elems[1] = v1; \
  this->elems[2] = v2; \
  this->elems[3] = v3; \
  this->elems[4] = v4; \
  this->elems[5] = v5; \
  this->elems[6] = v6; \
} \
class_name( \
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3, \
  const ElementType& v4, \
  const ElementType& v5, \
  const ElementType& v6, \
  const ElementType& v7 \
) : m_size(8) { \
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
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3, \
  const ElementType& v4, \
  const ElementType& v5, \
  const ElementType& v6, \
  const ElementType& v7, \
  const ElementType& v8 \
) : m_size(9) { \
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
  const ElementType& v0, \
  const ElementType& v1, \
  const ElementType& v2, \
  const ElementType& v3, \
  const ElementType& v4, \
  const ElementType& v5, \
  const ElementType& v6, \
  const ElementType& v7, \
  const ElementType& v8, \
  const ElementType& v9 \
) : m_size(10) { \
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
}
