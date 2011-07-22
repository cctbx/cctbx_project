#ifndef CCTBX_ADPTBX_HIRSHFELD_H
#define CCTBX_ADPTBX_HIRSHFELD_H

#include <cctbx/uctbx.h>
#include <scitbx/vec3.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/matrix/matrix_vector_operations.h>
#include <boost/array.hpp>

namespace cctbx { namespace adptbx {

  /// Mean square displacement along a direction
  /** Given a vector z and a displacement tensor U,
      both in fractional coordinates, the said quantity reads
      \f[
          H = \frac{ (Gz)^T U (Gz) }{ z^T G z }
      \f]
      where G is the metric matrix.

      A Taylor expansion reveals that the differential of H as a function
      of (G, z, U) reads
      \f[
           H'(G,z,U).(\gamma, \zeta, \upsilon) = \frac{1}{z^T G z} \big(
           (Gz)^T \upsilon (Gz)
           + 2 (GUGz - HGz)^T \zeta
           \\ + (2 UGz - Hz)^T \gamma z \big)
      \f]

      This class favours speed, the trade-off being a higher memory footprint
      and its keeping references to external objects.
   */
  template <typename T>
  class mean_square_displacement
  {
  public:
    /// Prepare a mean square displacement along the given direction
    mean_square_displacement(uctbx::unit_cell const &unit_cell,
                             scitbx::vec3<T> const &z)
    : unit_cell(unit_cell),
      g(unit_cell.metrical_matrix()),
      z(z),
      gz(g*z),
      d(gz*z)
    {
      if (!well_defined()) return;
      // compute the differential of H wrt u
      h_u[0] =   gz[0]*gz[0]; // H'_{U_11}
      h_u[1] =   gz[1]*gz[1]; // H'_{U_22}
      h_u[2] =   gz[2]*gz[2]; // H'_{U_33}
      h_u[3] = 2*gz[0]*gz[1]; // H'_{U_12}
      h_u[4] = 2*gz[0]*gz[2]; // H'_{U_13}
      h_u[5] = 2*gz[1]*gz[2]; // H'_{U_23}
      h_u /= d;
    }

    /// Compute the linearisation of the mean square displacement of u along z
    mean_square_displacement &operator()(scitbx::sym_mat3<T> const &u) {
      using scitbx::matrix::matrix_transposed_vector;

      if (!well_defined()) return *this;

      // compute H
      scitbx::vec3<T>ugz = u*gz;
      h = ugz*gz/d;


      h_z = 2.*(g*ugz - h*gz)/d;

      scitbx::vec3<T> tau = 2.*ugz - h*z;
      h_g[0] = tau[0]*z[0];               // H'_{g_11}
      h_g[1] = tau[1]*z[1];               // H'_{g_22}
      h_g[2] = tau[2]*z[2];               // H'_{g_33}
      h_g[3] = tau[0]*z[1] + tau[1]*z[0]; // H'_{g_12}
      h_g[4] = tau[0]*z[2] + tau[2]*z[0]; // H'_{g_13}
      h_g[5] = tau[1]*z[2] + tau[2]*z[1]; // H'_{g_23}
      h_g /= d;

      // H'_g G'_{a,b,c,alpha,beta,gamma}, i.e. [row] x [symmetric matrix]
      /* ugly but till someone writes a proper tiny matrix framework,
         as good as it gets
       */
      matrix_transposed_vector(6,6,
                               unit_cell.d_metrical_matrix_d_params().begin(),
                               h_g.begin(),
                               h_unit_cell_params.begin());
      return *this;
    }

    /// Whether the vector z is non-zero
    bool well_defined() { return d; }

    /// Value H
    T value() { return h;   }

    /// Gradient wrt U
    /** The scitbx::sym_mat3 ordering of matrix coefficients is assumed
     */
    af::tiny<T,6> const &grad_u() { return h_u; }

    /// Gradient wrt z
    af::tiny<T,3> const &grad_z() { return h_z; }

    /// Gradient wrt G
    /** The scitbx::sym_mat3 ordering of matrix coefficients is assumed
     */
    af::tiny<T,6> const &grad_g() { return h_g; }

    /// Gradient wrt unit cell parameters
    /** Usual order: \f$a, b, c, \alpha, \beta, \gamma\f$
     */
    af::tiny<T,6> const &grad_unit_cell_params() { return h_unit_cell_params; }

  private:
    uctbx::unit_cell const &unit_cell;
    scitbx::sym_mat3<T> const &g;
    scitbx::vec3<T> z;
    scitbx::vec3<T> gz;
    T d;
    T h;
    af::tiny<T,6> h_u, h_g, h_unit_cell_params;
    af::tiny<T,3> h_z;
  };

}}

#endif // Guard
