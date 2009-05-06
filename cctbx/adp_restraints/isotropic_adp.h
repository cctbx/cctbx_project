#ifndef CCTBX_ADP_RESTRAINTS_ISOTROPIC_ADP_H
#define CCTBX_ADP_RESTRAINTS_ISOTROPIC_ADP_H

#include <cctbx/error.h>
#include <cctbx/adptbx.h>
#include <scitbx/array_family/versa.h>
#include <vector>
#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace adp_restraints {

  struct isotropic_adp_proxy
  {
    //! Default constructor. Some data members are not initialized!
    isotropic_adp_proxy() {}

    //! Constructor.
    isotropic_adp_proxy(
      unsigned const& i_seq_,
      double weight_)
    :
      i_seq(i_seq_),
      weight(weight_)
    {}

    //! Indices into array of sites.
    unsigned i_seq;
    //! weight
    double weight;
  };

  class isotropic_adp
  {
  public:
    //! Default constructor. Some data members are not initialized!
    isotropic_adp() {}

    //! Constructor.
    isotropic_adp(
      scitbx::sym_mat3<double> const& u_cart_,
      double weight_)
    :
      u_cart(u_cart_),
      weight(weight_)
    {
      init_deltas();
    }

    //! Constructor.
    isotropic_adp(
      af::const_ref<scitbx::sym_mat3<double> > const& u_cart_,
      isotropic_adp_proxy const& proxy);

    af::tiny<double, 6>
    deltas() { return deltas_; }

    //! weight * [[sum_{ii} (deltas)**2] + [2 * sum_{i<j} (deltas)**2]].
    double residual() const;

    //! sqrt(mean_sq(deltas))
    //! The off-diagonal elements are included twice.
    double rms_deltas() const;

    //! This returns gradients_u_star
    scitbx::sym_mat3<double>
    gradients() const;

    //! Support for isotropic_adp_residual_sum.
    /*! Not available in Python.
     */
    void
    add_gradients(
      af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart,
      unsigned const& i_seqs) const;

    //! Cartesian anisotropic displacement parameters.
    scitbx::sym_mat3<double> u_cart;
    double weight;

  protected:
    void init_deltas();
    af::tiny<double, 6> deltas_;

  };

  /*! Fast computation of sum of isotropic_adp::residual() and gradients
      given an array of isotropic_adp proxies.
   */
  /*! The isotropic_adp::gradients() are added to the gradient_array if
      gradient_array.size() == sites_cart.size().
      gradient_array must be initialized before this function
      is called.
      No gradient calculations are performed if gradient_array.size() == 0.
   */
  double
  isotropic_adp_residual_sum(
    af::const_ref<scitbx::sym_mat3<double> > const& u_cart,
    af::const_ref<isotropic_adp_proxy> const& proxies,
    af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart);

  /*! \brief Fast computation of isotropic_adp::residual() given an array
      of adp_similarity proxies.
   */
  af::shared<double>
  isotropic_adp_residuals(
    af::const_ref<scitbx::sym_mat3<double> > const& u_cart,
    af::const_ref<isotropic_adp_proxy> const& proxies);

  /*! \brief Fast computation of isotropic_adp::rms_deltas() given an array
      of isotropic_adp proxies.
   */
  af::shared<double>
  isotropic_adp_deltas_rms(
    af::const_ref<scitbx::sym_mat3<double> > const& u_cart,
    af::const_ref<isotropic_adp_proxy> const& proxies);

}} // cctbx::adp_restraints

#endif
