#ifndef CCTBX_ADP_RESTRAINTS_ADP_SIMILARITY_H
#define CCTBX_ADP_RESTRAINTS_ADP_SIMILARITY_H

#include <cctbx/error.h>
#include <cctbx/adptbx.h>

#include <scitbx/array_family/versa.h>
#include <vector>
#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace adp_restraints {

  struct adp_similarity_proxy
  {
    //! Default constructor. Some data members are not initialized!
    adp_similarity_proxy() {}

    //! Constructor.
    adp_similarity_proxy(
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

  class adp_similarity
  {
  public:
    //! Default constructor. Some data members are not initialized!
    adp_similarity() {}

    //! Constructor.
    adp_similarity(
      af::tiny<scitbx::sym_mat3<double>, 2> const& u_cart_,
      af::tiny<double, 2> const& u_iso_,
      af::tiny<bool, 2> const& use_u_aniso_,
      double weight_)
    :
      u_cart(u_cart_),
      u_iso(u_iso_),
      use_u_aniso(use_u_aniso_),
      weight(weight_)
    {
      init_deltas();
    }

    //! Constructor.
    adp_similarity(
      af::const_ref<scitbx::sym_mat3<double> > const& u_cart_,
      af::const_ref<double> const& u_iso_,
      af::const_ref<bool> const& use_u_aniso_,
      adp_similarity_proxy const& proxy);

    af::tiny<double, 6>
    deltas() { return deltas_; }

    //! weight * sum(deltas)**2.
    double residual() const;

    scitbx::sym_mat3<double>
    gradient_0() const;

    //! This returns gradients_u_cart and gradients_u_equiv combined
    af::tiny<scitbx::sym_mat3<double>, 2>
    gradients() const;

    //! Support for adp_similarity_residual_sum.
    /*! Not available in Python.
     */
    void
    add_gradients(
      af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart,
      af::ref<double> const& gradients_iso,
      af::tiny<unsigned, 2> const& i_seqs) const;

    //! Cartesian anisotropic displacement parameters.
    af::tiny<scitbx::sym_mat3<double>, 2> u_cart;
    //! Isotropic displacement parameters.
    af::tiny<double, 2> u_iso;
    //! Use U aniso.
    af::tiny<bool, 2> use_u_aniso;

    double weight;

  protected:
    void init_deltas();
    af::tiny<double, 6> deltas_;

  };

  /*! Fast computation of sum of adp_similarity::residual() and gradients
      given an array of adp_similarity proxies.
   */
  /*! The adp_similarity::gradients() are added to the gradient_array if
      gradient_array.size() == sites_cart.size().
      gradient_array must be initialized before this function
      is called.
      No gradient calculations are performed if gradient_array.size() == 0.
   */
  double
  adp_similarity_residual_sum(
    af::const_ref<scitbx::sym_mat3<double> > const& u_cart,
    af::const_ref<double> const& u_iso,
    af::const_ref<bool> const& use_u_aniso,
    af::const_ref<adp_similarity_proxy> const& proxies,
    af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart,
    af::ref<double> const& gradients_iso);

}} // cctbx::adp_restraints

#endif
