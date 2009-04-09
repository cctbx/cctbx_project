#ifndef SCITBX_MATRIX_SPECIAL_MATRICES_H
#define SCITBX_MATRIX_SPECIAL_MATRICES_H

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/versa_matrix.h>

namespace scitbx { namespace matrix {

  template <typename NumType>
  af::versa<NumType, af::c_grid<2> > identity(int n) {
    af::versa<NumType, af::c_grid<2> > result(af::c_grid<2>(n,n));
    af::matrix_diagonal_set_in_place(result.ref(), NumType(1));
    return result;
  }

  template <typename NumType>
  af::versa<NumType, af::c_grid<2> > diagonal(af::const_ref<NumType> d) {
    std::size_t n = d.size();
    af::versa<NumType, af::c_grid<2> > result(af::c_grid<2>(n,n));
    af::ref<NumType, af::c_grid<2> > r = result.ref();
    for (int i=0; i<n; ++i) r(i,i) = d[i];
    return result;
  }

  template <typename NumType>
  af::versa<NumType, af::c_grid<2> > bidiagonal(
    af::const_ref<NumType> diagonal,
    af::const_ref<NumType> superdiagonal)
  {
    std::size_t n = diagonal.size();
    af::versa<NumType, af::c_grid<2> > result(af::c_grid<2>(n,n));
    af::ref<NumType, af::c_grid<2> > r = result.ref();
    for (int i=0; i<n; ++i) {
      r(i,i) = diagonal[i];
      if (i < n-1) r(i,i+1) = superdiagonal[i];
    }
    return result;
  }

}} // scitbx::matrix

#endif // GUARD
