#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_TAGS_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_TAGS_H

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

#define TRAITS(klass, periodic, padded) \
template <>                             \
struct may_be_padded<klass>             \
{                                       \
  static bool const value = periodic;   \
};                                      \
template <>                             \
struct may_be_periodic<klass>           \
{                                       \
  static bool const value = padded;     \
};

TRAITS(c_grid_padded<3>, true, false)
TRAITS(c_grid_periodic<3>, false, true);
TRAITS(c_grid_padded_periodic<3>, true, true)
TRAITS(flex_grid<>, true, false)

#undef TRAITS

}} // scitbx::af

#endif // GUARD
