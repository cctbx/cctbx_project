#ifndef SCITBX_MAT_REF_ROW_AND_COLUMN_H
#define SCITBX_MAT_REF_ROW_AND_COLUMN_H

#include <scitbx/mat_ref.h>
#include <scitbx/array_family/accessors/striding.h>

namespace scitbx {

template <typename ValueType>
af::ref<ValueType, af::striding_linear_accessor>
column_below(mat_ref<ValueType> const &a, std::size_t i, std::size_t j) {
  af::striding_linear_accessor acc(a.n_rows() - i, a.n_columns());
  return af::ref<ValueType, af::striding_linear_accessor>(&a(i,j), acc);
}

template <typename ValueType>
af::ref<ValueType>
row_right_of(mat_ref<ValueType> const &a, std::size_t i, std::size_t j) {
  return af::ref<ValueType>(&a(i,j), a.n_columns() - j);
}

} // scitbx

#endif // GUARD
