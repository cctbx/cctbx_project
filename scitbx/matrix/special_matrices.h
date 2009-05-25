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

  /// m x n matrix whose diagonal is d
  /** m and n shall not be smaller than the size p of d: if any of them is
      greater, then the first p elements of the diagonal are taken from d
      whereas the trailing part of the diagonal is zero
  */
  template <typename NumType>
  af::versa<NumType, af::c_grid<2> > diagonal(af::const_ref<NumType> d,
                                              int m, int n)
  {
    std::size_t p = d.size();
    SCITBX_ASSERT(m >= p)(m)(p);
    SCITBX_ASSERT(n >= p)(m)(p);
    af::versa<NumType, af::mat_grid> result(af::mat_grid(m,n));
    af::ref<NumType, af::mat_grid> r = result.ref();
    for (int i=0; i<p; ++i) r(i,i) = d[i];
    return result;
  }

  /// Overloaded variant returning a square diagonal matrix
  template <typename NumType>
  af::versa<NumType, af::c_grid<2> > diagonal(af::const_ref<NumType> d) {
    return diagonal(d, d.size(), d.size());
  }

  template <typename NumType>
  af::versa<NumType, af::c_grid<2> > upper_bidiagonal(
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

  template <typename NumType>
  af::versa<NumType, af::c_grid<2> > lower_bidiagonal(
    af::const_ref<NumType> diagonal,
    af::const_ref<NumType> subdiagonal)
  {
    std::size_t n = diagonal.size();
    af::versa<NumType, af::c_grid<2> > result(af::c_grid<2>(n,n));
    af::ref<NumType, af::c_grid<2> > r = result.ref();
    for (int i=0; i<n; ++i) {
      r(i,i) = diagonal[i];
      if (i < n-1) r(i+1,i) = subdiagonal[i];
    }
    return result;
  }

}} // scitbx::matrix

#endif // GUARD
