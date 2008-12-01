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
    unsigned n_rows;
    af::tiny<unsigned, NCols> col_perm;
    unsigned rank;
    unsigned nullity;

    full_pivoting_small() {}

    full_pivoting_small(
      af::ref<NumType, af::c_grid<2> > const& a_work,
      NumType const& min_abs_pivot=0,
      unsigned max_rank=NCols)
    {
      SCITBX_ASSERT(a_work.accessor()[0] <= MaxNRows);
      SCITBX_ASSERT(a_work.accessor()[1] == NCols);
      n_rows = a_work.accessor()[0];
      rank = full_pivoting_impl::reduction(
        n_rows,
        NCols,
        a_work.begin(),
        static_cast<NumType*>(0), // b
        min_abs_pivot,
        max_rank,
        col_perm.begin());
      nullity = NCols - rank;
    }

    af::tiny<NumType, NCols>
    back_substitution(
      const NumType* a_work,
      af::small<NumType, NCols> const& free_values) const
    {
      SCITBX_ASSERT(free_values.size() == nullity);
      af::tiny<NumType, NCols> perm_result;
      af::tiny<NumType, NCols> result;
      bool have_solution = full_pivoting_impl::back_substitution(
        n_rows,
        NCols,
        a_work,
        static_cast<const NumType*>(0),
        col_perm.begin(),
        rank,
        free_values.begin(),
        static_cast<NumType>(0), // epsilon
        perm_result.begin(),
        result.begin());
      SCITBX_ASSERT(have_solution);
      return result;
    }
  };

}}} // namespace scitbx::matrix::row_echelon

#endif // include guard
