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
    af::versa<NumType, af::flex_grid<> > a_work;
    af::shared<NumType> b_work;
    unsigned n_rows;
    af::shared<unsigned> col_perm;
    unsigned rank;
    unsigned nullity;

    full_pivoting() {}

    full_pivoting(
      af::versa<NumType, af::flex_grid<> > a_work_,
      NumType const& min_abs_pivot=0,
      int max_rank=-1)
    :
      a_work(a_work_)
    {
      reduction_(min_abs_pivot, max_rank);
    }

    full_pivoting(
      af::versa<NumType, af::flex_grid<> > a_work_,
      af::shared<NumType> b_work_,
      NumType const& min_abs_pivot=0,
      int max_rank=-1)
    :
      a_work(a_work_),
      b_work(b_work_)
    {
      reduction_(min_abs_pivot, max_rank);
    }

    protected:
    void
    reduction_(
      NumType const& min_abs_pivot,
      int max_rank)
    {
      if (a_work.accessor().nd() != 2) {
        throw std::runtime_error("a_work matrix must be two-dimensional.");
      }
      af::c_grid<2> c_grid(a_work.accessor());
      n_rows = static_cast<unsigned>(c_grid[0]);
      unsigned n_cols = static_cast<unsigned>(c_grid[1]);
      col_perm.resize(n_cols);
      rank = full_pivoting_impl::reduction(
        n_rows,
        n_cols,
        a_work.begin(),
        (b_work.size() == 0 ? static_cast<NumType*>(0) : b_work.begin()),
        min_abs_pivot,
        (max_rank < 0 ? n_cols : static_cast<unsigned>(max_rank)),
        col_perm.begin());
      nullity = n_cols - rank;
    }
    public:

    bool
    is_in_row_space(
      af::const_ref<NumType> const& x,
      NumType const& epsilon) const
    {
      SCITBX_ASSERT(x.size() == col_perm.size());
      af::shared<NumType> x_copy(x.begin(), x.end());
      return full_pivoting_impl::is_in_row_space(
        static_cast<unsigned>(col_perm.size()),
        a_work.begin(),
        col_perm.begin(),
        rank,
        x_copy.begin(),
        epsilon);
    }

    boost::optional<af::shared<NumType> >
    back_substitution(
      af::const_ref<NumType> const& free_values,
      NumType const& epsilon=0) const
    {
      SCITBX_ASSERT(free_values.size() == nullity);
      af::shared<NumType> perm_result(col_perm.size());
      af::shared<NumType> result(col_perm.size());
      bool have_solution = full_pivoting_impl::back_substitution(
        n_rows,
        static_cast<unsigned>(col_perm.size()),
        a_work.begin(),
        (b_work.size() == 0 ? static_cast<const NumType*>(0) : b_work.begin()),
        col_perm.begin(),
        rank,
        free_values.begin(),
        epsilon,
        perm_result.begin(),
        result.begin());
      if (!have_solution) return boost::optional<af::shared<NumType> >();
      return boost::optional<af::shared<NumType> >(result);
    }
  };

}}} // namespace scitbx::matrix::row_echelon

#endif // include guard
