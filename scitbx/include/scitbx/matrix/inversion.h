#ifndef SCITBX_MATRIX_INVERSION_H
#define SCITBX_MATRIX_INVERSION_H

#include <boost/scoped_array.hpp>
#include <algorithm>
#include <stdexcept>

namespace scitbx { namespace matrix {

  //! Generic in-place matrix inversion with pivoting.
  /*! On input a is a pointer to the values of a matrix with n rows and
      n columns. b is a pointer to the values of m column vectors with
      n elements per vector.
      On output a is the inverted matrix and b contains the solutions
      x of a*x=b.
      A std::runtime_error is thrown if the input matrix is singular.
   */
  template <typename FloatType>
  void
  inversion_in_place(
    FloatType *a,
    std::size_t n,
    FloatType *b,
    std::size_t m)
  {
    if (n == 0) return;
    const std::size_t max_n_stack = 10;
    std::size_t scratch_stack[max_n_stack*3];
    boost::scoped_array<std::size_t> scratch_dynamic;
    std::size_t* ipivot;
    if (n <= max_n_stack) {
      ipivot = scratch_stack;
    }
    else {
      boost::scoped_array<std::size_t> scratch(new std::size_t[n*3]);
      scratch_dynamic.swap(scratch);
      ipivot = scratch_dynamic.get();
    }
    std::fill(ipivot, ipivot+n, 0);
    std::size_t* indxc = ipivot + n;
    std::size_t* indxr = indxc + n;
    for(std::size_t i=0;i<n;i++) {
      FloatType a_abs_max = 0;
      std::size_t irow = 0;
      std::size_t icol = 0;
      for(std::size_t j=0;j<n;j++) {
        if (ipivot[j] != 1) {
          for(std::size_t k=0;k<n;k++) {
            if (ipivot[k] == 0) {
              FloatType ajk_abs = a[j*n+k];
              if (ajk_abs < 0) ajk_abs = -ajk_abs;
              if (ajk_abs >= a_abs_max) {
                a_abs_max = ajk_abs;
                irow = j;
                icol = k;
              }
            }
            else if (ipivot[k] > 1) {
              throw std::runtime_error("inversion_in_place: singular matrix");
            }
          }
        }
      }
      ipivot[icol]++;
      if (irow != icol) {
        for(std::size_t l=0;l<n;l++) std::swap(a[irow*n+l], a[icol*n+l]);
        for(std::size_t l=0;l<m;l++) std::swap(b[l*n+irow], b[l*n+icol]);
      }
      indxr[i] = irow;
      indxc[i] = icol;
      if (a[icol*n+icol] == 0) {
        throw std::runtime_error("inversion_in_place: singular matrix");
      }
      FloatType one_over_pivot = 1 / a[icol*n+icol];
      a[icol*n+icol] = 1;
      for(std::size_t l = 0; l < n; l++) a[icol*n+l] *= one_over_pivot;
      for(std::size_t l = 0; l < m; l++) b[l*n+icol] *= one_over_pivot;
      for(std::size_t ll = 0; ll < n; ll++) {
        if (ll != icol) {
          FloatType f = a[ll*n+icol];
          a[ll*n+icol] = 0;
          for(std::size_t l=0;l<n;l++) a[ll*n+l] -= a[icol*n+l] * f;
          for(std::size_t l=0;l<m;l++) b[l*n+ll] -= b[l*n+icol] * f;
        }
      }
    }
    for(std::size_t l=n;l>0;) {
      l--;
      if (indxr[l] != indxc[l]) {
        for(std::size_t k=0;k<n;k++) {
          std::swap(a[k*n+indxr[l]], a[k*n+indxc[l]]);
        }
      }
    }
  }

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_INVERSION_H
