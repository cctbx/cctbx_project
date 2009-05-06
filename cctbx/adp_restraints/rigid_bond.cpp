#include <assert.h>
#include <math.h>
#include <iostream>
#include <cctbx/adptbx.h>
#include <cctbx/error.h>
#include <cctbx/adp_restraints/rigid_bond.h>

using namespace std;
namespace cctbx { namespace adp_restraints {

  rigid_bond_pair::rigid_bond_pair(vec3<double> const& site1,
                                   vec3<double> const& site2,
                                   sym_mat3<double> const& ustar1,
                                   sym_mat3<double> const& ustar2,
                                   cctbx::uctbx::unit_cell const& uc)
  /* Hirshfeld's Rigid Bond Test (1: Acta Cryst. (1976). A32, 239,
                                  2: SHELX manual)
  */
  {
    sym_mat3<double> g = uc.metrical_matrix();
    vec3<double> l_12 = site1 - site2;
    vec3<double> l_21 = site2 - site1;
    double bond_length_sq = l_12 * g * l_12;
    z_12_ = (g * l_12) * ustar1 * (g * l_12) / bond_length_sq;
    z_21_ = (g * l_21) * ustar2 * (g * l_21) / bond_length_sq;
    delta_z_ = std::abs(z_12_ - z_21_);
  }

  //! Constructor.
  rigid_bond::rigid_bond(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<scitbx::sym_mat3<double> > const& u_cart_,
    rigid_bond_proxy const& proxy)
  :
    weight(proxy.weight)
  {
    CCTBX_ASSERT(sites_cart.size() == u_cart_.size());
    for (int i=0;i<2;i++) {
      std::size_t i_seq = proxy.i_seqs[i];
      CCTBX_ASSERT(i_seq < sites_cart.size());
      sites[i] = sites_cart[i_seq];
      u_cart[i] = u_cart_[i_seq];
    }
    init_delta();
  }

  void
  rigid_bond::init_delta()
  {
    l_12 = sites[0] - sites[1];
    vec3<double> l_21 = -l_12;
    bond_length_sq = scitbx::fn::pow2(l_12.length());
    z_12_ = l_12 * u_cart[0] * l_12 / bond_length_sq;
    z_21_ = l_21 * u_cart[1] * l_21 / bond_length_sq;
    delta_z_ = std::abs(z_12_ - z_21_);
  }

  //! Gradient of residual with respect to u_cart[0]
  scitbx::sym_mat3<double>
  rigid_bond::gradient_0() const
  {
    scitbx::sym_mat3<double> result;
    for (int i=0;i<3;i++) {
      result[i] = scitbx::fn::pow2(l_12[i]);
    }
    result[3] = 2 * l_12[0] * l_12[1];
    result[4] = 2 * l_12[0] * l_12[2];
    result[5] = 2 * l_12[1] * l_12[2];
    result *= -2 * weight * delta_z_ / bond_length_sq;
    return result;
  }

  af::tiny<scitbx::sym_mat3<double>, 2>
  rigid_bond::gradients() const
  {
    af::tiny<scitbx::sym_mat3<double>, 2> result;
    result[0] = gradient_0();
    result[1] = -result[0];
    return result;
  }

  void
  rigid_bond::add_gradients(
    af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart,
    af::tiny<unsigned, 2> const& i_seqs) const
  {
    scitbx::sym_mat3<double> g0 = gradient_0();
    gradients_aniso_cart[i_seqs[0]] += g0;
    gradients_aniso_cart[i_seqs[1]] += -g0;
  }

  double
  rigid_bond_residual_sum(
  af::const_ref<scitbx::vec3<double> > const& sites_cart,
  af::const_ref<scitbx::sym_mat3<double> > const& u_cart,
  af::const_ref<rigid_bond_proxy> const& proxies,
  af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart)
  {
    CCTBX_ASSERT(   gradients_aniso_cart.size() == 0
                 || gradients_aniso_cart.size() == sites_cart.size());
    double result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      rigid_bond_proxy const& proxy = proxies[i];
      rigid_bond restraint(sites_cart, u_cart, proxy);
      result += restraint.residual();
      if (gradients_aniso_cart.size() != 0) {
        restraint.add_gradients(gradients_aniso_cart, proxy.i_seqs);
      }
    }
    return result;
  }

  af::shared<double>
  rigid_bond_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<scitbx::sym_mat3<double> > const& u_cart,
    af::const_ref<rigid_bond_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      rigid_bond_proxy const& proxy = proxies[i];
      rigid_bond restraint(sites_cart, u_cart, proxy);
      result.push_back(restraint.residual());
    }
    return result;
  }

  af::shared<double>
  rigid_bond_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<scitbx::sym_mat3<double> > const& u_cart,
    af::const_ref<rigid_bond_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      rigid_bond_proxy const& proxy = proxies[i];
      rigid_bond restraint(sites_cart, u_cart, proxy);
      result.push_back(restraint.delta_z());
    }
    return result;
  }

}} // namespace cctbx::apd_restraints
