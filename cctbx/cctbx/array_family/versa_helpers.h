// $Id$
// Included from versa_plain.h
// DO NOT INCLUDE THIS FILE DIRECTLY.

#define CCTBX_ARRAY_FAMILY_VERSA_CONSTRUCTORS(class_name) \
class_name() {} \
class_name(const accessor_type& ac, const ElementType& x = ElementType()) { \
  this->resize(ac); \
  std::fill(this->begin(), this->end(), x); \
} \
explicit \
class_name(int n0) { \
  this->resize(accessor_type(n0)); \
  std::fill(this->begin(), this->end(), ElementType()); \
} \
class_name(int n0, int n1) { \
  this->resize(accessor_type(n0, n1)); \
  std::fill(this->begin(), this->end(), ElementType()); \
} \
class_name(int n0, int n1, int n2) { \
  this->resize(accessor_type(n0, n1, n2)); \
  std::fill(this->begin(), this->end(), ElementType()); \
}
