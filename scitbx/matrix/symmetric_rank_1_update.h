#ifndef SCITBX_MATRIX_SYMMETRIC_RANK_1_UPDATE_H
#define SCITBX_MATRIX_SYMMETRIC_RANK_1_UPDATE_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <scitbx/matrix/matrix_vector_operations.h>

namespace scitbx { namespace matrix {

  template <typename T>
  class sum_of_symmetric_rank_1_updates
  {
  private:
    af::versa<T, af::packed_u_accessor> a;

  public:
    sum_of_symmetric_rank_1_updates(int n)
    : a(n)
    {}

    void add(af::const_ref<double> const &x, double alpha) {
      SCITBX_ASSERT(x.size() == a.accessor().n_rows())(x.size())(a.accessor().n_rows());
      add(x.begin(), alpha);
    }

    void add(double const *x, double alpha) {
      symmetric_packed_u_rank_1_update(a.accessor().n_rows(), a.begin(), x, alpha);
    }

    void reset() {
      std::fill(a.begin(), a.end(), T(0));
    }

    void finalise() {}

    operator af::versa<double, af::packed_u_accessor>() {
      return a;
    }
  };
}}


#endif
