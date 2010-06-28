#ifndef CCTBX_ADP_RESTRAINTS_RIGID_BOND_H
#define CCTBX_ADP_RESTRAINTS_RIGID_BOND_H

#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <cctbx/adptbx.h>

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

  class rigid_bond
  {
  public:
    //! Default constructor. Some data members are not initialized!
    rigid_bond() {}

    //! Constructor.
    rigid_bond(
      af::tiny<scitbx::vec3<double>, 2> const& sites_,
      af::tiny<scitbx::sym_mat3<double>, 2> const& u_cart_,
      double weight_)
    :
      sites(sites_),
      u_cart(u_cart_),
      weight(weight_)
    {
      init_delta();
    }

    //! Constructor.
    rigid_bond(
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

    //! weight * delta_z**2.
    double
    residual() const { return weight * scitbx::fn::pow2(delta_z_); }

    //! Gradient of delta_z with respect to u_cart[0]
    scitbx::sym_mat3<double>
    grad_delta_0() const
    {
      scitbx::sym_mat3<double> result;
      for (int i=0;i<3;i++) {
        result[i] = scitbx::fn::pow2(l_12[i]);
      }
      result[3] = 2 * l_12[0] * l_12[1];
      result[4] = 2 * l_12[0] * l_12[2];
      result[5] = 2 * l_12[1] * l_12[2];
      result *= -1 / bond_length_sq;
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

    double z_12() { return z_12_; }
    double z_21() { return z_21_; }
    double delta_z() { return delta_z_; }

    //! Cartesian coordinates of bonded sites.
    af::tiny<scitbx::vec3<double>, 2> sites;
    //! Cartesian anisotropic displacement parameters.
    af::tiny<scitbx::sym_mat3<double>, 2> u_cart;

    double delta;
    double weight;
  protected:
    void init_delta()
    {
      l_12 = sites[0] - sites[1];
      vec3<double> l_21 = -l_12;
      bond_length_sq = scitbx::fn::pow2(l_12.length());
      z_12_ = l_12 * u_cart[0] * l_12 / bond_length_sq;
      z_21_ = l_21 * u_cart[1] * l_21 / bond_length_sq;
      delta_z_ = std::abs(z_12_ - z_21_);
    }

    double z_12_, z_21_, delta_z_;
    //! The vector in the direction of the bond site 1 -> site 2
    vec3<double> l_12;
    double bond_length_sq;
  };

  /*! Fast computation of sum of rigid_bond::residual() and gradients
      given an array of rigid_bond proxies.
   */
  /*! The rigid_bond::gradients() are added to the gradient_array if
      gradient_array.size() == sites_cart.size().
      gradient_array must be initialized before this function
      is called.
      No gradient calculations are performed if gradient_array.size() == 0.
   */
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

  /*! \brief Fast computation of rigid_bond::residual() given an array
      of rigid_bond proxies.
   */
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

  /*! \brief Fast computation of rigid_bond::delta_z() given an array
      of rigid_bond proxies.
   */
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

}} // namespace cctbx::adp_restraints

#endif // CCTBX_ADP_RESTRAINTS_RIGID_BOND_H
