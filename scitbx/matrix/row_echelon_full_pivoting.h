#ifndef SCITBX_MATRIX_ROW_ECHELON_FULL_PIVOTING_H
#define SCITBX_MATRIX_ROW_ECHELON_FULL_PIVOTING_H

#include <scitbx/matrix/row_echelon_full_pivoting_impl.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace scitbx { namespace matrix { namespace row_echelon {

  template <typename NumType>
  struct full_pivoting
  {
    af::versa<NumType, af::c_grid<2> > echelon_form;
    af::shared<unsigned> row_perm;
    af::shared<unsigned> col_perm;
    unsigned n_pivots;

    full_pivoting() {}

    full_pivoting(
      af::versa<NumType, af::flex_grid<> >& matrix,
      NumType const& min_abs_pivot=0,
      int max_rank=-1)
    {
      af::c_grid<2> c_grid(matrix.accessor());
      unsigned n_rows = static_cast<unsigned>(c_grid[0]);
      unsigned n_cols = static_cast<unsigned>(c_grid[1]);
      row_perm.resize(n_rows);
      col_perm.resize(n_cols);
      n_pivots = full_pivoting_impl::reduction(
        matrix.begin(),
        n_rows,
        n_cols,
        min_abs_pivot,
        (max_rank < 0 ? n_cols : static_cast<unsigned>(max_rank)),
        row_perm.begin(),
        col_perm.begin());
      c_grid = af::c_grid<2>(n_pivots, n_cols);
      matrix.resize(c_grid.as_flex_grid());
      echelon_form = af::versa<double, af::c_grid<2> >(
        matrix.handle(), c_grid);
    }

    unsigned
    rank() const { return n_pivots; }

    unsigned
    nullity() const
    {
      return static_cast<unsigned>(col_perm.size()) - n_pivots;
    }

    bool
    is_in_row_span(
      af::const_ref<NumType> const& vector,
      NumType const& epsilon) const
    {
      SCITBX_ASSERT(vector.size() == col_perm.size());
      af::shared<NumType> vector_copy(vector.begin(), vector.end());
      return full_pivoting_impl::is_in_row_span(
        static_cast<unsigned>(col_perm.size()),
        echelon_form.begin(),
        col_perm.begin(),
        n_pivots,
        vector_copy.begin(),
        epsilon);
    }

    af::shared<NumType>
    back_substitution(af::const_ref<NumType> const& free_values) const
    {
      SCITBX_ASSERT(free_values.size() == nullity());
      af::shared<NumType> perm_result(col_perm.size());
      af::shared<NumType> result(col_perm.size());
      full_pivoting_impl::back_substitution(
        static_cast<unsigned>(col_perm.size()),
        echelon_form.begin(),
        col_perm.begin(),
        n_pivots,
        free_values.begin(),
        perm_result.begin(),
        result.begin());
      return result;
    }
  };

}}} // namespace scitbx::matrix::row_echelon

#endif // include guard
