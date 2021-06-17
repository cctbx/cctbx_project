#ifndef SCITBX_ARRAY_FAMILY_FLEX_TYPES_H
#define SCITBX_ARRAY_FAMILY_FLEX_TYPES_H

#include <complex>
#include <stdint.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>

namespace scitbx { namespace af {

  template <typename ElementType>
  struct flex
  {
    typedef versa<ElementType, flex_grid<> > type;
  };

  template <typename ElementType>
  struct flex_const_ref
  {
    typedef const_ref<ElementType, flex_grid<> > type;
  };

  template <typename ElementType>
  struct flex_ref
  {
    typedef ref<ElementType, flex_grid<> > type;
  };

#define SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(T, N) \
  typedef flex<T >::type flex_ ## N; \
  typedef flex_const_ref<T >::type flex_ ## N ## _const_ref; \
  typedef flex_ref<T >::type flex_ ## N ## _ref;

  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(bool, bool)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(int, int)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(long, long)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(std::size_t, size_t)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(float, float)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(double, double)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(std::complex<double>, complex_double)

  // fixed width integer types
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(int8_t, int8_t)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(int16_t, int16_t)
  // SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(int32_t, int32_t)
  #if defined(_MSC_VER)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(int64_t, int64_t)
  #endif
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(uint8_t, uint8_t)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(uint16_t, uint16_t)
  SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(uint32_t, uint32_t)
  // SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS(uint64_t, uint64_t)

#undef SCITBX_ARRAY_FAMILY_FLEX_TYPEDEFS

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_FLEX_TYPES_H
