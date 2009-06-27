#ifndef SCITBX_MATRIX_NORMS_H
#define SCITBX_MATRIX_NORMS_H

#include <scitbx/array_family/accessors/mat_grid.h>
#include <cmath>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/math/accumulators.h>

namespace scitbx { namespace matrix {

template <typename FloatType>
FloatType norm_1(af::const_ref<FloatType, af::mat_grid> const &a) {
  af::shared<FloatType> sum_col_(a.n_columns(), 0);
  af::ref<FloatType> sum_col = sum_col_.ref();
  for (int i=0; i<a.n_rows(); ++i) for (int j=0; j<a.n_columns(); ++j) {
    sum_col[j] += std::abs(a(i,j));
  }
  return af::max(sum_col);
}

template <typename FloatType>
FloatType norm_1(af::const_ref<FloatType> const &x) {
  FloatType result = 0;
  for (int i=0; i<x.size(); ++i) result += std::abs(x[i]);
  return result;
}

template <typename FloatType>
FloatType norm_inf(af::const_ref<FloatType, af::mat_grid> const &a) {
  FloatType result = 0;
  for (int i=0; i<a.n_rows(); ++i) {
    FloatType sum_row = 0;
    for (int j=0; j<a.n_columns(); ++j) {
      sum_row += std::abs(a(i,j));
    }
    result = std::max(result, sum_row);
  }
  return result;
}

template <typename FloatType>
FloatType norm_frobenius(af::const_ref<FloatType, af::mat_grid> const &a) {
  af::const_ref<FloatType> a_ij = a.as_1d();
  math::accumulator::norm_accumulator<FloatType> norm_accu;
  for (int k=0; k < a_ij.size(); ++k) norm_accu(a_ij[k]);
  return norm_accu.norm();
}

}} // scitbx::matrix

#endif // GUARD
