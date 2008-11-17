#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_TRAITS_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_TRAITS_H

#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/accessors/c_grid_periodic.h>
#include <scitbx/array_family/accessors/c_grid_padded_periodic.h>
#include <scitbx/array_family/accessors/flex_grid.h>

namespace scitbx { namespace af {

template <class GridType>
struct may_be_padded
{
  static bool const value = false;
};

template <class GridType>
struct may_be_periodic
{
  static bool const value = false;
};

#define SCITBX_LOC(klass, periodic, padded) \
template <>                             \
struct may_be_padded<klass >            \
{                                       \
  static bool const value = periodic;   \
};                                      \
template <>                             \
struct may_be_periodic<klass >          \
{                                       \
  static bool const value = padded;     \
};

SCITBX_LOC(c_grid_padded<3>, true, false)
SCITBX_LOC(c_grid_periodic<3>, false, true)
SCITBX_LOC(c_grid_padded_periodic<3>, true, true)
SCITBX_LOC(flex_grid<>, true, false)

#undef SCITBX_LOC

}} // scitbx::af

#endif // GUARD
