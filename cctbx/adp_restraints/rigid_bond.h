#ifndef CCTBX_ADP_RESTRAINTS_RIGID_BOND_H
#define CCTBX_ADP_RESTRAINTS_RIGID_BOND_H

#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <cctbx/adptbx.h>
#include <cctbx/adp_restraints/adp_restraints.h>
#include <scitbx/matrix/matrix_vector_operations.h>

namespace cctbx { namespace adp_restraints {

using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;

  /* Hirshfeld's Rigid Bond Test (1: Acta Cryst. (1976). A32, 239,
                                  2: SHELX manual)
  */
  class rigid_bond_pair {
  public:
    rigid_bond_pair(vec3<double> const& site1,
                    vec3<double> const& site2,
                    sym_mat3<double> const& ustar1,
                    sym_mat3<double> const& ustar2,
                    cctbx::uctbx::unit_cell const& uc)
  {
    sym_mat3<double> g = uc.metrical_matrix();
    vec3<double> l_12 = site1 - site2;
    vec3<double> l_21 = site2 - site1;
    double bond_length_sq = l_12 * g * l_12;
    z_12_ = (g * l_12) * ustar1 * (g * l_12) / bond_length_sq;
    z_21_ = (g * l_21) * ustar2 * (g * l_21) / bond_length_sq;
    delta_z_ = std::abs(z_12_ - z_21_);
  }

    double z_12() { return z_12_; }
    double z_21() { return z_21_; }
    double delta_z() { return delta_z_; }

  private:
      double z_12_, z_21_, delta_z_;
  };

  struct rigid_bond_proxy
  {
    //! Default constructor. Some data members are not initialized!
    rigid_bond_proxy() {}

    //! Constructor.
    rigid_bond_proxy(
      af::tiny<unsigned, 2> const& i_seqs_,
      double weight_)
    :
      i_seqs(i_seqs_),
      weight(weight_)
    {}

    //! Indices into array of sites.
    af::tiny<unsigned, 2> i_seqs;
    //! weight
    double weight;
  };

  class rigid_bond {
  public:
    //! Constructor.
    rigid_bond(
      af::tiny<scitbx::vec3<double>, 2> const& sites,
      af::tiny<scitbx::sym_mat3<double>, 2> const& u_cart,
      double weight_)
    :
      weight(weight_)
    {
      init_delta(sites, u_cart);
    }

    //! Constructor.
    rigid_bond(
      adp_restraint_params<double> const &params,
      rigid_bond_proxy const& proxy)
    :
      weight(proxy.weight)
    {
      CCTBX_ASSERT(params.sites_cart.size() == params.u_cart.size());
      CCTBX_ASSERT(proxy.i_seqs[0] < params.sites_cart.size());
      CCTBX_ASSERT(proxy.i_seqs[1] < params.sites_cart.size());
      init_delta(
        af::tiny<scitbx::vec3<double>, 2>(
          params.sites_cart[proxy.i_seqs[0]], params.sites_cart[proxy.i_seqs[1]]),
        af::tiny<scitbx::sym_mat3<double>, 2>(
          params.u_cart[proxy.i_seqs[0]], params.u_cart[proxy.i_seqs[1]]));
    }

    //! weight * delta_z**2.
    double
    residual() const { return weight * scitbx::fn::pow2(delta_z_); }

    //! Gradient of delta_z with respect to u_cart[0]
    scitbx::sym_mat3<double> grad_delta_0() const {
      scitbx::sym_mat3<double> result;
      for (int i=0;i<3;i++) {
        result[i] = scitbx::fn::pow2(l_12[i]);
      }
      result[3] = 2 * l_12[0] * l_12[1];
      result[4] = 2 * l_12[0] * l_12[2];
      result[5] = 2 * l_12[1] * l_12[2];
      result /= bond_length_sq;
      return result;
    }

    //! Gradient of residual with respect to u_cart[0]
    scitbx::sym_mat3<double>
    gradient_0() const
    {
      scitbx::sym_mat3<double> result = grad_delta_0();
      result *= 2 * weight * delta_z_;
      return result;
    }

    af::tiny<scitbx::sym_mat3<double>, 2>
    gradients() const
    {
      af::tiny<scitbx::sym_mat3<double>, 2> result;
      result[0] = gradient_0();
      result[1] = -result[0];
      return result;
    }

    //! Support for rigid_bond_residual_sum.
    /*! Not available in Python.
     */
    void
    add_gradients(
      af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart,
      af::tiny<unsigned, 2> const& i_seqs) const
    {
      scitbx::sym_mat3<double> g0 = gradient_0();
      gradients_aniso_cart[i_seqs[0]] += g0;
      gradients_aniso_cart[i_seqs[1]] += -g0;
    }

    void
    linearise(
      uctbx::unit_cell const &unit_cell,
      cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
      cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
      af::tiny<unsigned, 2> const& i_seqs) const
    {
      af::const_ref<double, af::mat_grid> const &f
        = unit_cell.u_star_to_u_cart_linear_map();
      scitbx::sym_mat3<double> grad_u_cart = grad_delta_0();
      scitbx::sym_mat3<double> grad_u_star;
      scitbx::matrix::matrix_transposed_vector(
        6, 6, f.begin(), grad_u_cart.begin(), grad_u_star.begin());
      std::size_t row_i = linearised_eqns.next_row();
      for (std::size_t i=0;i<2;i++) {
        if (i == 1) grad_u_star = -grad_u_star;
        cctbx::xray::parameter_indices const &ids_i
          = parameter_map[i_seqs[i]];
        if (ids_i.u_aniso == -1) continue;
        for (std::size_t j=0;j<6;j++) {
          linearised_eqns.design_matrix(row_i, ids_i.u_aniso+j)
            = grad_u_star[j];
        }
      linearised_eqns.weights[row_i] = weight;
      linearised_eqns.deltas[row_i] = delta_z_;
      }
    }

    double z_12() { return z_12_; }
    double z_21() { return z_21_; }
    double delta_z() { return delta_z_; }

    double weight;
  protected:
    void init_delta(af::tiny<scitbx::vec3<double>, 2> const &sites,
      af::tiny<scitbx::sym_mat3<double>, 2> const &u_cart)
    {
      l_12 = sites[0] - sites[1];
      vec3<double> l_21 = -l_12;
      bond_length_sq = scitbx::fn::pow2(l_12.length());
      z_12_ = l_12 * u_cart[0] * l_12 / bond_length_sq;
      z_21_ = l_21 * u_cart[1] * l_21 / bond_length_sq;
      delta_z_ = z_12_ - z_21_;
    }

    double z_12_, z_21_, delta_z_;
    //! The vector in the direction of the bond site 1 -> site 2
    vec3<double> l_12;
    double bond_length_sq;
  };

  /*! \brief Fast computation of rigid_bond::delta_z() given an array
      of rigid_bond proxies.
   */
  af::shared<double>
  rigid_bond_deltas(
    adp_restraint_params<double> const &params,
    af::const_ref<rigid_bond_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0; i<proxies.size(); i++) {
      result.push_back(rigid_bond(params, proxies[i]).delta_z());
    }
    return result;
  }

}} // namespace cctbx::adp_restraints

#endif // CCTBX_ADP_RESTRAINTS_RIGID_BOND_H
