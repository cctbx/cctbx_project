#ifndef SCITBX_MATRIX_SYMMETRIC_RANK_1_UPDATE_H
#define SCITBX_MATRIX_SYMMETRIC_RANK_1_UPDATE_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <scitbx/matrix/matrix_vector_operations.h>
#include <scitbx/matrix/vector_operations.h>
#include <scitbx/array_family/simple_io.h>
#include <fast_linalg/lapacke.h>

namespace scitbx { namespace matrix {

  /// Sum of symmetric rank-1 updates \f$\alpha_i x_i x_i^T\f$
  template <typename T>
  class sum_of_symmetric_rank_1_updates
  {
  private:
    af::versa<T, af::packed_u_accessor> a;

  public:
    /// Initialise the sum to a zero matrix of size n
    sum_of_symmetric_rank_1_updates(int n)
    : a(n)
    {}

    /// Add \f$\alpha x x^T\f$ to the sum
    void add(af::const_ref<T> const &x, T alpha) {
      SCITBX_ASSERT(x.size() == a.accessor().n_rows())
      (x.size())
      (a.accessor().n_rows());
      add(x.begin(), alpha);
    }

    /// Add \f$\alpha x x^T\f$ (overload without size checks for speed)
    void add(T const *x, T alpha) {
      symmetric_packed_u_rank_1_update(a.accessor().n_rows(), a.begin(), x, alpha);
    }

    /// Add s in-place
    sum_of_symmetric_rank_1_updates
    &operator+=(sum_of_symmetric_rank_1_updates const &s) {
      af::ref<T> a_r = a.ref().as_1d();
      a_r += s.a.const_ref().as_1d();
      return *this;
    }

    /// Cancel all the rank-1 updates
    /** The sum is reset to the zero matrix */
    void reset() {
      std::fill(a.begin(), a.end(), T(0));
    }

    /// Called after after all rank-1 updates have been performed
    /** This does nothing in this class */
    void finalise() {}

    /// The resulting (symmetric) matrix
    /** It returns a meaningful result only after finalise() has been called */
    operator af::versa<T, af::packed_u_accessor>() {
      return a;
    }
  };

  /// Symmetric rank-N update \f$A^T A\f$, specified row by row
  /// but computed at BLAS 3 speed
  /** Thus this class is equivalent to class sum_of_symmetric_rank_1_updates
   *  when all \f$\alpha_i\f$ are non-negative, since then
   *  \f$\alpha x x^T = y y^T\f$ where \f$y = \sqrt{\alpha} x\f$ is one row
   *  of matrix A.
   */
  template <typename T>
  class rank_n_update
  {
  public:
    /// Prepare for a resulting matrix of size n
    rank_n_update(int n)
    : a((af::reserve(n*n/2))), aaT_rfp(n), aaT_packed(n), cols(n)
    {}

    /// Add a row \f$\sqrt{\alpha} x\f$ to matrix A
    /// Precondition: alpha >= 0
    void add(af::const_ref<T> const &x, T alpha) {
      SCITBX_ASSERT(x.size() == cols)(x.size())(cols);
      add(x.begin(), alpha);
    }

    /// Overload without size check for speed but alpha >= 0 still enforced
    void add(T const *x, T alpha) {
      SCITBX_ASSERT(alpha >= 0)(alpha);
      a.extend(x, x + cols);
      matrix::scale_vector(cols, a.end()-cols, std::sqrt(alpha));
    }

    /// Add u in-place
    rank_n_update &operator+=(rank_n_update const &u) {
      a.extend(u.a.begin(), u.a.end());
      return *this;
    }

    /// Cancel all the rank-1 updates
    void reset() {
      a.clear();
    }

    /// Called after after all rank-1 updates have been performed
    void finalise() {
      std::size_t const rows = a.size()/cols;
      // A^T A = [A^T] [A^T]^T and A^T is the column-major version of
      // the row-major A
      using namespace fast_linalg;
      sfrk(LAPACK_COL_MAJOR, 'N', 'L', 'N',
           cols, rows, 1.0, a.begin(), cols, 0.0, aaT_rfp.begin());
      int info = tfttp(LAPACK_COL_MAJOR, 'N', 'L',
                       cols, aaT_rfp.begin(), aaT_packed.begin());
      SCITBX_ASSERT(info == 0)(info);
    }

    /// The resulting (symmetric) matrix
    /** It returns a meaningful result only after finalise() has been called */
    operator af::versa<T, af::packed_u_accessor>() {
      return aaT_packed;
    }

  private:
    af::shared<T> a;
    af::versa<T, af::packed_u_accessor> aaT_rfp, aaT_packed;
    int cols;
  };

}}


#endif
