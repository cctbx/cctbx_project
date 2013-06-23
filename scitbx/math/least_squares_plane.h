#ifndef SCITBX_MATH_LEAST_SQUARE_PLANE_H
#define SCITBX_MATH_LEAST_SQUARE_PLANE_H

#include <scitbx/vec3.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/math/accumulators.h>
#include <scitbx/matrix/eigensystem.h>

namespace scitbx { namespace math {

/** "Best" plane through a set of points in dimension 3.

The plane H implementd by this class is such that H minimises
\f$ \sum_i d_i^2 \f$
where \f$ d_i \f$ is the distance of the i-th point to H. The

There is only one solution which is related to the first and
second statistical moments of the set of points: the barycentre \f$ b \f$
and the covariance tensor \f$ C \f$. Denoting by \f$ p_i \f$ the
column-vector of coordinates of the i-th point, they are defined as

- \f$ b = \frac{1}{n} \sum_i p_i \f$
- \f$ C = \frac{1}{n} \sum_i (p_i - b)(p_i - b)^T \f$

where \f$ n \f$ is the number of points.

Then H is the hyperplane passing through \f$ b \f$ and whose normal \f$ n \f$ is
the eigenvector corresponding to the smallest eigenvalue of \f$ C \f$.
Thus the equation of the hyperplane is \f$ n^T p = d \f$
where \f$ d = n^T b \f$ is the distance of H to the origin (this is
the so-called Hessian normal form of the equation of H).
*/
template<typename ScalarType=double>
class least_squares_plane
{
  public:
    typedef ScalarType scalar_type;
    typedef vec3<scalar_type> vector_type;
    typedef mat3<scalar_type> matrix_type;

    /// Construct for a set of points with the given 3-vector of coordinates
    least_squares_plane(af::const_ref<vector_type> const &points) {
      scitbx::math::accumulator::inertia_accumulator<scalar_type> acc;
      for(int i=0; i<points.size(); ++i) acc(points[i]);
      eigen_decomposition_type eigen(acc.covariance_matrix());
      // eigenvalues are sorted in decreasing order
      scalar_type *v = eigen.vectors().begin();
      n = vector_type(v[6], v[7], v[8]).normalize();
      d = n*acc.center_of_mass();
    }

    /// Normal \f$ n \f$ to the plane
    vector_type const &normal() { return n; }

    /// Distance \f$ d \f$ of the origin to the plane
    scalar_type distance_to_origin() { return d; }

  private:
    typedef scitbx::matrix::eigensystem::real_symmetric<scalar_type>
            eigen_decomposition_type;
    vector_type n;
    scalar_type d;
};

}}

#endif // GUARD
