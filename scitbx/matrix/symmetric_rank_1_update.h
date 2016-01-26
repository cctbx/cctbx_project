#ifndef SCITBX_MATRIX_SYMMETRIC_RANK_1_UPDATE_H
#define SCITBX_MATRIX_SYMMETRIC_RANK_1_UPDATE_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <scitbx/matrix/matrix_vector_operations.h>

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


}}


#endif
