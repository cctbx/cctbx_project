#ifndef SCITBX_MATRIX_LU_DECOMPOSITION_H
#define SCITBX_MATRIX_LU_DECOMPOSITION_H

#include <boost/scoped_array.hpp>
#include <algorithm>
#include <stdexcept>

namespace scitbx { namespace matrix {

  template <typename FloatType>
  void
  lu_decomposition_in_place(
    FloatType *a,
    std::size_t n,
    std::size_t *pivot_indices)
  {
    const std::size_t max_n_stack = 16;
    FloatType scratch_stack[max_n_stack];
    boost::scoped_array<FloatType> scratch_dynamic;
    FloatType *vv;
    if (n <= max_n_stack) {
      vv = scratch_stack;
    }
    else {
      boost::scoped_array<FloatType> scratch(new FloatType[n]);
      scratch_dynamic.swap(scratch);
      vv = scratch_dynamic.get();
    }
    pivot_indices[n] = 0;
    for(std::size_t i=0;i<n;i++) {
      FloatType big = 0;
      for(std::size_t j=0;j<n;j++) {
        FloatType dum = a[i*n+j];
        if (dum < 0) dum = -dum;
        if (dum > big) big = dum;
      }
      if (big == 0) {
        throw std::runtime_error("lu_decomposition_in_place: singular matrix");
      }
      vv[i] = 1 / big;
    }
    std::size_t imax = 0;
    for(std::size_t j=0;j<n;j++) {
      for(std::size_t i=0;i<j;i++) {
        FloatType sum = a[i*n+j];
        for(std::size_t k=0;k<i;k++) sum -= a[i*n+k] * a[k*n+j];
        a[i*n+j] = sum;
      }
      FloatType big = 0;
      for(std::size_t i=j;i<n;i++) {
        FloatType sum = a[i*n+j];
        for(std::size_t k=0;k<j;k++) sum -= a[i*n+k] * a[k*n+j];
        a[i*n+j] = sum;
        if (sum < 0) sum = -sum;
        FloatType dum = vv[i] * sum;
        if (dum >= big) {
          big = dum;
          imax = i;
        }
      }
      if (j != imax) {
        for(std::size_t k=0;k<n;k++) std::swap(a[imax*n+k], a[j*n+k]);
        pivot_indices[n]++;
        vv[imax] = vv[j]; /* no swap, we don't need vv[j] any more */
      }
      pivot_indices[j] = imax;
      if (a[j*n+j] == 0) {
        throw std::runtime_error("lu_decomposition_in_place: singular matrix");
      }
      if (j+1 < n) {
        FloatType dum = 1 / a[j*n+j];
        for(std::size_t i=j+1;i<n;i++)
          a[i*n+j] *= dum;
      }
    }
  }

  template <typename FloatType>
  void
  lu_back_substitution(
    const FloatType *a,
    std::size_t n,
    const std::size_t *pivot_indices,
    FloatType *b)
  {
    std::size_t ii = n;
    for(std::size_t i=0;i<n;i++) {
      std::size_t pivot_indices_i = pivot_indices[i];
      if (pivot_indices_i >= n) {
        throw std::runtime_error(
          "lu_back_substitution: pivot_indices[i] out of range");
      }
      FloatType sum = b[pivot_indices_i];
      b[pivot_indices_i] = b[i];
      if (ii != n) {
        for(std::size_t j=ii;j+1<=i;j++) sum -= a[i*n+j] * b[j];
      }
      else if (sum) {
        ii = i;
      }
      b[i] = sum;
    }
    for(std::size_t i=n;i>0;) {
      i--;
      FloatType sum = b[i];
      for(std::size_t j=i+1;j<n;j++) sum -= a[i*n+j] * b[j];
      b[i] = sum / a[i*n+i];
    }
  }

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_LU_DECOMPOSITION_H
