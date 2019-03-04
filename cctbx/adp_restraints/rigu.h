#ifndef CCTBX_ADP_RESTRAINTS_RIGU_H
#define CCTBX_ADP_RESTRAINTS_RIGU_H

#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <cctbx/adptbx.h>
#include <cctbx/adp_restraints/adp_restraints.h>
#include <scitbx/matrix/matrix_vector_operations.h>

namespace cctbx { namespace adp_restraints {

  /** RIGU restraint
   *  Restrains U33, U13 and U23 of 2 adps expressed in a cartesian base
   *  along the bond of the 2 atoms. U33 is aligned along the bond
   *  See Thorn, A., et. al. (2012). Acta Cryst. A68, 448-451.
   *  and Parois, P., et. al. (2018). J. Appl. Cryst. 51, 1059-1068.
   */

  struct rigu_proxy {
    //! Default constructor. Some data members are not initialized!
    rigu_proxy() {}

    //! Constructor.
    rigu_proxy(
      af::tiny<unsigned, 2> const& i_seqs_, double weight_)
    : i_seqs(i_seqs_),
      weight(weight_)
    {}

    //! Indices into array of sites.
    af::tiny<unsigned, 2> i_seqs;
    //! weight
    double weight;
  };

  class rigu {
  public:
    typedef scitbx::mat3<double> mat3d;
    typedef scitbx::sym_mat3<double> sym_mat3d;
    typedef scitbx::vec3<double> vec3d;

    //! Constructor.
    rigu(
      af::tiny<vec3d, 2> const& sites,
      af::tiny<sym_mat3d, 2> const& u_cart,
      double weight_)
    :
      dRUcart(3),
      weight(weight_)
    {
      init_delta(sites, u_cart);
      calc_gradients();
    }

    //! Constructor.
    rigu(
      adp_restraint_params<double> const &params,
      rigu_proxy const& proxy)
    :
      dRUcart(3),
      weight(proxy.weight)
    {
      CCTBX_ASSERT(params.sites_cart.size() == params.u_cart.size());
      CCTBX_ASSERT(proxy.i_seqs[0] < params.sites_cart.size());
      CCTBX_ASSERT(proxy.i_seqs[1] < params.sites_cart.size());
      init_delta(
        af::tiny<vec3d, 2>(
          params.sites_cart[proxy.i_seqs[0]],
          params.sites_cart[proxy.i_seqs[1]]),
        af::tiny<sym_mat3d, 2>(
          params.u_cart[proxy.i_seqs[0]],
          params.u_cart[proxy.i_seqs[1]]));
      calc_gradients();
    }

    //! weight * delta[i]**2.
    double
    residual33() const {
      return weight * scitbx::fn::pow2(delta_33_);
    }
    double
    residual13() const {
      return weight * scitbx::fn::pow2(delta_13_);
    }
    double
    residual23() const {
      return weight * scitbx::fn::pow2(delta_23_);
    }
    double
    residual() const {
      return residual33() + residual13() + residual23();
    }

    //! Gradient of residual with respect to u_cart[0]
    sym_mat3d gradient_33() const {
      return dRUcart[0] * (2 * weight * delta_33_);
    }
    sym_mat3d gradient_13() const {
      return dRUcart[1] * (2 * weight * delta_13_);
    }
    sym_mat3d gradient_23() const {
      return dRUcart[2] * (2 * weight * delta_23_);
    }

    af::tiny<sym_mat3d, 2> gradients33() const {
      af::tiny<scitbx::sym_mat3<double>, 2> result;
      result[0] = gradient_33();
      result[1] = -result[0];
      return result;
    }

    af::tiny<sym_mat3d, 2> gradients13() const {
      af::tiny<scitbx::sym_mat3<double>, 2> result;
      result[0] = gradient_13();
      result[1] = -result[0];
      return result;
    }

    af::tiny<sym_mat3d, 2> gradients23() const {
      af::tiny<scitbx::sym_mat3<double>, 2> result;
      result[0] = gradient_23();
      result[1] = -result[0];
      return result;
    }

    //! Support for rigu_residual_sum.
    /*! Not available in Python.
     */
    void add_gradients(
      af::ref<sym_mat3d> const& gradients_aniso_cart,
      af::tiny<unsigned, 2> const& i_seqs) const
    {
      sym_mat3d g0 = gradient_33();
      gradients_aniso_cart[i_seqs[0]] += g0;
      gradients_aniso_cart[i_seqs[1]] += -g0;
      g0 = gradient_13();
      gradients_aniso_cart[i_seqs[0]] += g0;
      gradients_aniso_cart[i_seqs[1]] += -g0;
      g0 = gradient_23();
      gradients_aniso_cart[i_seqs[0]] += g0;
      gradients_aniso_cart[i_seqs[1]] += -g0;
    }

    void linearise(
      uctbx::unit_cell const &unit_cell,
      cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
      cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
      af::tiny<unsigned, 2> const& i_seqs) const
    {
      af::const_ref<double, af::mat_grid> const &f
        = unit_cell.u_star_to_u_cart_linear_map();
      scitbx::sym_mat3<double> grad_u_star;
      double const deltas[3] = {delta_33_, delta_13_, delta_23_};

      for (std::size_t k = 0; k < 3; k++) {
        scitbx::matrix::matrix_transposed_vector(
          6, 6, f.begin(), dRUcart[k].begin(), grad_u_star.begin());
        std::size_t row_i = linearised_eqns.next_row();
        for (std::size_t i = 0; i < 2; i++) {
          if (i == 1) {
            grad_u_star = -grad_u_star;
          }
          cctbx::xray::parameter_indices const &ids_i
            = parameter_map[i_seqs[i]];
          if (ids_i.u_aniso == -1) {
            continue;
          }
          for (std::size_t j = 0; j < 6; j++) {
            linearised_eqns.design_matrix(row_i, ids_i.u_aniso+j)
              = grad_u_star[j] * (j < 3 ? 1 : 2);
          }
          linearised_eqns.weights[row_i] = weight;
          linearised_eqns.deltas[row_i] = deltas[k];
        }
      }
    }

    // for testing
    mat3d getRM() const {
      return RM;
    }
    // for testing
    af::shared<sym_mat3d> raw_gradients() const {
      return dRUcart;
    }

    double delta_33() const { return delta_33_; }
    double delta_13() const { return delta_13_; }
    double delta_23() const { return delta_23_; }
    double delta() const { return delta_33_+delta_13_+delta_23_; }

  protected:
    void init_delta(af::tiny<vec3d, 2> const &sites,
      af::tiny<sym_mat3d, 2> const &u_cart)
    {

      /** Calculating the rotation matrix to align U33 along the bond
       *  z-axis along the bond
       */
      vec3d rot3 = sites[1] - sites[0];

      //! Any perpendicular axis will do
      vec3d rot2(rot3[2], rot3[2], -rot3[0]-rot3[1]);

      if(rot2.length_sq() < 1e-4) {
        rot2[0] = -rot3[1]-rot3[2];
        rot2[1] = rot3[1];
        rot2[2] = rot3[1];
      }

      //! Last axis to form a direct orthonormal basis
      vec3d rot1 = rot2.cross(rot3);

      RM.set_row(0, rot1.normalize());
      RM.set_row(1, rot2.normalize());
      RM.set_row(2, rot3.normalize());
      mat3d RMt = RM.transpose();

      //! Calculating Ucart in the new basis
      mat3d RUcart1 = (RM*u_cart[0])*RMt;
      mat3d RUcart2 = (RM*u_cart[1])*RMt;

      //! calculating the 3 deltas involved
      delta_33_ = RUcart1(2,2) - RUcart2(2,2);
      delta_13_ = RUcart1(0,2) - RUcart2(0,2);
      delta_23_ = RUcart1(1,2) - RUcart2(1,2);

      /** Update weight see
       *  Thorn, A., et. al. (2012). Acta Cryst. A68, 448-451
       *  sqrt(p**2 + UeqA + UeqB) * d * p *sigma
       */
      double d = rot3.length();
      double const p = 0.5;
      double UeqA = u_cart[0].trace()/3.0;
      double UeqB = u_cart[1].trace()/3.0;
      weight = weightfinal(weight, p, d, UeqA, UeqB);
    }

    void calc_gradients() {
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j <= i; j++) {
          dRUcart[0](i, j) = RM(2, i)*RM(2, j);
          dRUcart[1](i, j) = RM(0, i)*RM(2, j);
          dRUcart[2](i, j) = RM(1, i)*RM(2, j);
          if (i != j) {
            dRUcart[0](i, j) += RM(2, j)*RM(2, i);
            dRUcart[1](i, j) += RM(0, j)*RM(2, i);
            dRUcart[2](i, j) += RM(1, j)*RM(2, i);
          }
        }
      }
    }

    /** weights in RIGU depends on the distance between atoms,
     *  Ueqs and a fixed parameter p
     */
    double weightfinal(double weight, double p, double d,
      double UeqA, double UeqB) const
    {
      return p*p / ((p*p + UeqA + UeqB) * d*d) * weight;
    }

    double delta_33_;
    double delta_13_;
    double delta_23_;
    mat3d RM;
    af::shared<sym_mat3d> dRUcart;
  public:
    double weight;
  };

}} // namespace cctbx::adp_restraints

#endif // CCTBX_ADP_RESTRAINTS_RIGU_H
