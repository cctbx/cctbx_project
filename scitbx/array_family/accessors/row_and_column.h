#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_ROW_AND_COLUMN_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_ROW_AND_COLUMN_H

#include <scitbx/array_family/accessors/mat_grid.h>
#include <scitbx/array_family/accessors/striding.h>

namespace scitbx { namespace af {

template <typename ValueType>
af::ref<ValueType, af::striding_linear_accessor>
column_below(af::ref<ValueType, af::mat_grid> const &a,
             std::size_t i, std::size_t j)
{
  af::striding_linear_accessor acc(a.n_rows() - i, a.n_columns());
  return af::ref<ValueType, af::striding_linear_accessor>(&a(i,j), acc);
}

template <typename ValueType>
af::ref<ValueType>
row_right_of(af::ref<ValueType, af::mat_grid> const &a,
             std::size_t i, std::size_t j)
{
  return af::ref<ValueType>(&a(i,j), a.n_columns() - j);
}

template <typename ValueType>
af::const_ref<ValueType>
row(af::const_ref<ValueType, af::mat_grid> const &a, std::size_t i) {
  return af::const_ref<ValueType>(&a(i, 0), a.n_columns());
}

}} // scitbx::af

#endif // GUARD
