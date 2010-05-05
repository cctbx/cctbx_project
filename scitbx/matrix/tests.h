#ifndef SCITBX_MATRIX_TESTS_H
#define SCITBX_MATRIX_TESTS_H

#include <scitbx/array_family/accessors/mat_grid.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/array_family/ref_algebra.h>
#include <scitbx/matrix/norms.h>
#include <scitbx/matrix/packed.h>
#include <limits>
#include <cmath>

namespace scitbx { namespace matrix {

  /// A measure of whether the columns and/or the rows of u
  /// form an orthonormal system.
  /** Reference: DLQT02.F and DQRT02.F in LAPACK tests (c.f. very end of those
      files) */
  template <typename T>
  T normality_ratio(af::const_ref<T, af::mat_grid> const &u,
                    T eps=std::numeric_limits<T>::epsilon())
  {
    typedef af::versa<T, af::mat_grid> matrix_t;
    typedef af::const_ref<T, af::mat_grid> matrix_const_ref_t;
    typedef af::ref<T, af::mat_grid> matrix_ref_t;

    int m=u.n_rows(), n = u.n_columns();
    matrix_t ut_= af::matrix_transpose(u);
    matrix_const_ref_t ut = ut_.const_ref();
    if (m <= n) { // only the rows can be tested
      matrix_t delta = af::matrix_multiply(u, ut);
      matrix_ref_t delta_ = delta.ref();
      for (int i=0; i<m; ++i) delta_(i,i) -= 1;
      return (norm_1(delta.const_ref())/n)/eps;
    }
    else { // only the columns can be tested
      matrix_t delta = af::matrix_multiply(ut, u);
      matrix_ref_t delta_ = delta.ref();
      for (int i=0; i<n; ++i) delta_(i,i) -= 1;
      return (norm_1(delta.const_ref())/m)/eps;
    }
  }


  /// A measure the relative difference between a and b
  template <typename T>
  T equality_ratio(af::const_ref<T, af::mat_grid> const &a,
                   af::const_ref<T, af::mat_grid> const &b,
                   T eps=std::numeric_limits<T>::epsilon())
  {
    SCITBX_ASSERT(a.n_rows() == b.n_rows());
    SCITBX_ASSERT(a.n_columns() == b.n_columns());
    typedef af::c_grid<2> dim;
    int m=a.n_rows(), n=a.n_columns();
    af::versa<T, dim> delta(dim(m, n));
    for (int i=0; i<m; ++i) for (int j=0; j<n; ++j) delta(i,j) = a(i,j) - b(i,j);
    return ((norm_1(delta.const_ref())
             /std::max(a.n_rows(), a.n_columns()))
             /norm_1(a))/eps;
  }

  /// Test of the solution of A x = b using the Cholesky decomposition A = R^T R
  template <typename T>
  T cholesky_test_ratio(af::const_ref<T, af::mat_grid> const &a,
                        af::const_ref<T> const &x,
                        af::const_ref<T> const &b,
                        T eps=std::numeric_limits<T>::epsilon())
  {
    af::shared<T> y = matrix_multiply(a, x);
    af::shared<T> delta = y.ref() - b;
    return norm_1(delta.ref())/(norm_1(a) * norm_1(x) * eps);
  }

  /// Test the quality of the inverse of a obtained from A Cholesky decomposition
  /** Terribly inefficient implementation */
  template <typename T>
  T residual_of_symmetric(af::const_ref<T, af::packed_u_accessor> const &a,
                          af::const_ref<T, af::packed_u_accessor> const &a_inv,
                          T eps=std::numeric_limits<T>::epsilon())
  {
    SCITBX_ASSERT(a.accessor().n == a_inv.accessor().n);
    std::size_t n = a.accessor().n;
    af::versa<T, af::mat_grid> a_ = packed_u_as_symmetric(a.as_1d()),
                               a_inv_ = packed_u_as_symmetric(a_inv.as_1d());
    af::versa<T, af::mat_grid> delta = af::matrix_multiply(a_.ref(),
                                                           a_inv_.ref());
    af::matrix_diagonal_add_in_place(delta.ref(), T(-1));
    return norm_1(delta.ref())/(n*norm_1(a_.ref())*norm_1(a_inv_.ref())*eps);
  }


}} // scitbx::matrix

#endif // GUARD
