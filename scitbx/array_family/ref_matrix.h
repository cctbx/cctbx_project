#ifndef SCITBX_ARRAY_FAMILY_REF_MATRIX_H
#define SCITBX_ARRAY_FAMILY_REF_MATRIX_H

/**  Matrix related functions operating on af::ref and af::const_ref. */

#include <scitbx/array_family/ref.h>
#include <scitbx/matrix/multiply.h>

namespace scitbx { namespace af {

  //! Matrix multiplication.
  template <typename NumTypeA,  typename AccessorTypeA,
            typename NumTypeB,  typename AccessorTypeB,
            typename NumTypeAB, typename AccessorTypeAB>
  inline
  void
  multiply(
    const_ref<NumTypeA, AccessorTypeA> const& a,
    const_ref<NumTypeB, AccessorTypeB> const& b,
    ref<NumTypeAB, AccessorTypeAB> const& ab)
  {
    SCITBX_ASSERT(a.n_columns() == b.n_rows());
    SCITBX_ASSERT(ab.n_rows() == a.n_rows());
    SCITBX_ASSERT(ab.n_columns() == b.n_columns());
    matrix::multiply(a.begin(), b.begin(),
                     a.n_rows(), a.n_columns(), b.n_columns(),
                     ab.begin());
  }

}}

#endif // GUARD
