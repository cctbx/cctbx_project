#ifndef SCITBX_MATRIX_ROW_ECHELON_FULL_PIVOTING_SMALL_H
#define SCITBX_MATRIX_ROW_ECHELON_FULL_PIVOTING_SMALL_H

#include <scitbx/matrix/row_echelon_full_pivoting_impl.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace scitbx { namespace matrix { namespace row_echelon {

  template <typename NumType, unsigned MaxNRows, unsigned NCols>
  struct full_pivoting_small
  {
    af::small<NumType, NCols*NCols> echelon_form;
    af::small<unsigned, MaxNRows> row_perm;
    af::tiny<unsigned, NCols> col_perm;
    af::small<unsigned, NCols> pivot_cols;
    af::small<unsigned, NCols> free_cols;

    full_pivoting_small() {}

    full_pivoting_small(
      af::ref<NumType, af::c_grid<2> > const& matrix,
      NumType const& min_abs_pivot=0,
      unsigned max_rank=NCols)
    {
      SCITBX_ASSERT(matrix.accessor()[0] <= MaxNRows)(matrix.accessor()[0])
                                                     (MaxNRows);
      SCITBX_ASSERT(matrix.accessor()[1] == NCols);
      unsigned n_rows = matrix.accessor()[0];
      row_perm.resize(n_rows);
      pivot_cols.resize(NCols);
      free_cols.resize(NCols);
      unsigned pivot_cols_size = full_pivoting_impl::reduction(
        matrix.begin(),
        n_rows,
        NCols,
        min_abs_pivot,
        max_rank,
        row_perm.begin(),
        col_perm.begin(),
        pivot_cols.begin(),
        free_cols.begin());
      pivot_cols.resize(pivot_cols_size);
      free_cols.resize(NCols - pivot_cols_size);
      // copy result to local memory
      echelon_form.assign(matrix.begin(),
                          matrix.begin() + pivot_cols_size*NCols);
    }

    unsigned
    row_rank() const
    {
      return static_cast<unsigned>(pivot_cols.size());
    }

    bool
    is_in_row_span(
      af::small<NumType, NCols> vector, // copy
      NumType const& epsilon) const
    {
      SCITBX_ASSERT(vector.size() == NCols);
      return full_pivoting_impl::is_in_row_span(
        NCols,
        echelon_form.begin(),
        col_perm.begin(),
        pivot_cols.begin(),
        static_cast<unsigned>(pivot_cols.size()),
        vector.begin(),
        epsilon);
    }

    af::tiny<NumType, NCols>
    back_substitution(af::small<NumType, NCols> const& free_values) const
    {
      SCITBX_ASSERT(free_values.size() == free_cols.size());
      af::tiny<NumType, NCols> perm_result;
      af::tiny<NumType, NCols> result;
      full_pivoting_impl::back_substitution(
        NCols,
        echelon_form.begin(),
        col_perm.begin(),
        pivot_cols.begin(),
        static_cast<unsigned>(pivot_cols.size()),
        free_cols.begin(),
        static_cast<unsigned>(free_cols.size()),
        free_values.begin(),
        perm_result.begin(),
        result.begin());
      return result;
    }
  };

}}} // namespace scitbx::matrix::row_echelon

#endif // include guard
