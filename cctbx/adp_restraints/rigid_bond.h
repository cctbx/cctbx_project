#ifndef CCTBX_ADP_RESTRAINTS_RIGID_BOND_H
#define CCTBX_ADP_RESTRAINTS_RIGID_BOND_H

#include <cctbx/error.h>
#include <cctbx/adptbx.h>

#include <scitbx/array_family/versa.h>
#include <vector>
#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace adp_restraints {

using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;

class rigid_bond_pair {
public:
    rigid_bond_pair(vec3<double> const& site1,
                    vec3<double> const& site2,
                    sym_mat3<double> const& ustar1,
                    sym_mat3<double> const& ustar2,
                    cctbx::uctbx::unit_cell const& uc);
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
      rigid_bond_proxy const& proxy);

    //! weight * delta_z**2.
    double
    residual() const { return weight * scitbx::fn::pow2(delta_z_); }

    scitbx::sym_mat3<double>
    gradient_0() const;

    af::tiny<scitbx::sym_mat3<double>, 2>
    gradients() const;

    //! Support for rigid_bond_residual_sum.
    /*! Not available in Python.
     */
    void
    add_gradients(
      af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart,
      af::tiny<unsigned, 2> const& i_seqs) const;

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
    void init_delta();

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
    af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart);

}} // namespace cctbx::apd_restraints

#endif // CCTBX_ADP_RESTRAINTS_RIGID_BOND_H
