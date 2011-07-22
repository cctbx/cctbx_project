#ifndef CCTBX_ADPTBX_HIRSHFELD_H
#define CCTBX_ADPTBX_HIRSHFELD_H

#include <scitbx/vec3.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/matrix/matrix_vector_operations.h>
#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/rt_mx.h>

// Implementation
#include <scitbx/sparse/vector.h>

namespace cctbx { namespace adptbx {

  namespace details {
    template <typename T>
    class sparse_grad_container : public af::small<T, 2*(3+6)>
    {};
  }

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

  /// Relative difference of the mean square displacements of two scatterers
  /** Given two scatterers at sites \f$x_1, x_2\f$ with ADP \f$u_1, u_2\f$,
   the relative difference is defined as
   \f[ h = 2 \frac{h_1 - h_2}{h_1 + h_2} \f]
   where \f$h_1, h_2\f$ are the mean square displacement along \f$x_1 - x_2\f$
   of respectively \f$u_1, u_2\f$.
   */
  template <typename T>
  class relative_hirshfeld_difference {
    T val;
    af::tiny<T, 3> h_x_1, h_x_2;
    af::tiny<T,6> h_u_1, h_u_2, h_unit_cell_params;

  public:
    /// Construct for the given scatterer
    /** The 2nd site and displacement tensor may be the image through \f$r_2\f$
     of the given site and displacement \f$x_2, u_2\f$.
     It is necessary to know about this for correct a computation of the gradient.
     */
    relative_hirshfeld_difference(uctbx::unit_cell const &uc,
                                  scitbx::vec3<T> const &x1,
                                  scitbx::sym_mat3<T> const &u1,
                                  scitbx::vec3<T> const &x2,
                                  scitbx::sym_mat3<T> const &u2,
                                  sgtbx::rt_mx const &r_2)
    {
      using scitbx::matrix::matrix_transposed_vector;

      // Compute mean square displacement linearisation for scatterer 1 and 2
      scitbx::vec3<T> z = x1 - r_2(x2);
      adptbx::mean_square_displacement<T> msd(uc, z);
      msd(u1);
      T h1 = msd.value();
      af::tiny<T, 3> h1_z = msd.grad_z();
      af::tiny<T, 6> h1_u = msd.grad_u();
      af::tiny<T, 6> h1_uc = msd.grad_unit_cell_params();
      msd(r_2(u2));
      T h2 = msd.value();
      af::tiny<T, 3> const &h2_z = msd.grad_z();
      af::tiny<T, 6> const &h2_u = msd.grad_u();
      af::tiny<T, 6> const &h2_uc = msd.grad_unit_cell_params();

      // Compute linearisation of h
      T d = 1./(h1 + h2), d_sq = d*d;
      val = 2*(h1 - h2)*d;
      T h_h1 = 4*h2*d_sq, h_h2 = -4*h1*d_sq;

      // h'_{x_1} =  h'_z, h'_{x_2} = -h'_z r_2
      h_x_1 = h_h1*h1_z + h_h2*h2_z;
      h_x_2 = -h_x_1;
      if (!r_2.is_unit_mx()) {
        h_x_2 = h_x_2 * r_2.r().as_double();
      }

      /* h'_{u_1} = h'_{h_1} h'_{1; u}
         h'_{u_2} = h'_{h_2} h'_{2; u} R_2
         where R_2 is the 6x6 matrix transforming u_2 into r_2 u_2 r_2^T
         when those symmetric tensors are represented
         as coefficients (11,22,33,12,13,23)
       */
      h_u_1 = h_h1*h1_u;
      h_u_2 = h_h2*h2_u;
      if (!r_2.is_unit_mx()) {
        af::tiny<T, 6*6> rr_2 = r_2.r().tensor_transform_matrix<T>();
        af::tiny<T, 6> rhs;
        matrix_transposed_vector(6,6,
                                 rr_2.begin(), h_u_2.begin(), rhs.begin());

        h_u_2 = rhs;
      }

      // h'_{a,b,c,alpha,beta,gamma} = h'_{h_1} h'_{1;...} + h'_{h_2} h'_{2;...}
      h_unit_cell_params = h_h1*h1_uc + h_h2*h2_uc;
    }

    /// The value of h
    T value() { return val; }

    /// Differential of h wrt to \f$x_1\f$
    af::tiny<T, 3> const &grad_x1() const { return h_x_1; }

    /// Differential of h wrt to \f$x_2\f$
    af::tiny<T, 3> const &grad_x2() const { return h_x_2; }

    /// Differential of h wrt to \f$u_1\f$
    af::tiny<T, 6> const &grad_u1() const { return h_u_1; }

    /// Differential of h wrt to \f$u_2\f$
    af::tiny<T, 6> const &grad_u2() const { return h_u_2; }

    /// Differential of h wrt to unit cell parameters
    af::tiny<T, 6> const &grad_unit_cell_params() const {
      return h_unit_cell_params;
    }

    /// Estimated standard deviation of h
    /** Given the variance matrix for crystallographic parameters and
        the indices of relevant sites and ADP's therein
     */
    T esd(af::const_ref<T, af::packed_u_accessor> const &variance,
          std::size_t i_x1, std::size_t i_u1,
          std::size_t i_x2, std::size_t i_u2,
          af::tiny<T, 6> const &unit_cell_param_sigmas)
    {
      // sites and ADP's uncertainties
      scitbx::sparse::vector<T, details::sparse_grad_container>
      h_prime(variance.n_rows());
      for (int i=0; i<3; i++) h_prime[i_x1 + i] = h_x_1[i];
      for (int i=0; i<3; i++) h_prime[i_x2 + i] = h_x_2[i];
      for (int i=0; i<6; i++) h_prime[i_u1 + i] = h_u_1[i];
      for (int i=0; i<6; i++) h_prime[i_u2 + i] = h_u_2[i];
      T result = scitbx::sparse::quadratic_form(h_prime, variance, h_prime);

      // unit cell parameter uncertainties
      for (int i=0; i<6; i++) {
        result += h_unit_cell_params[i]*unit_cell_param_sigmas[i];
      }

      return std::sqrt(result);
    }

  };

  /*template <typename T>
   af::shared< af::tiny<T,2> >
   relative_differences(uctbx::unit_cell const &unit_cell,
   af::const_ref<scatterer<T> > const &scatterers,
   parameter_map<scatterer<T> > const &map,
   crystal::pair_sym_table const &pair_sym_table,
   af::const_ref<T, af::packed_u_accessor> const &variance)
   {
   CCTBX_ASSERT(scatterers.size() == map.size())
   (scatterers.size())(map.size());
   CCTBX_ASSERT(scatterers.size() == pair_sym_table.size())
   (scatterers.size())(pair_sym_table.size());
   for (std::size_t i=0; i<pair_sym_table.size(); i++) {
   scatterer<T> const &sc_i = scatterers[i];
   BOOST_FOREACH (crystal::pair_sym_dict::value_type const &idx_op,
   pair_sym_table[i])
   {
   std::size_t j = idx_op.first;
   scatterer<T> const &sc_j = scatterers[j];
   BOOST_FOREACH (crystal::pair_sym_ops::value_type const &r,
   idx_op.second) {
   adptbx::mean_square_displacement<T> msd(sc_i.site, r(sc_j.site));
   msd(sc_i.u_star);

   }
   }

   }
   }
   */


}}

#endif // Guard
