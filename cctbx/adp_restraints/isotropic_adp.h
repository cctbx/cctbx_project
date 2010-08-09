#ifndef CCTBX_ADP_RESTRAINTS_ISOTROPIC_ADP_H
#define CCTBX_ADP_RESTRAINTS_ISOTROPIC_ADP_H

#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <cctbx/adptbx.h>
#include <cctbx/restraints.h>

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
      isotropic_adp_proxy const& proxy)
    :
      weight(proxy.weight)
    {
      for (int i=0;i<2;i++) {
        std::size_t i_seq = proxy.i_seq;
        CCTBX_ASSERT(i_seq < u_cart_.size());
        u_cart = u_cart_[i_seq];
      }
      init_deltas();
    }


    scitbx::sym_mat3<double>
    deltas() { return deltas_; }

    //! weight * [[sum_{ii} (deltas)**2] + [2 * sum_{i<j} (deltas)**2]].
    /* This is the square of the Frobenius norm of the matrix of deltas, or
       alternatively the inner product of the matrix of deltas with itself.
       This is used since the residual is then rotationally invariant.
     */
    double residual() const
    {
      return weight * deltas_.dot(deltas_);
    }

    //! sqrt(mean_sq(deltas))
    //! The off-diagonal elements are included twice.
    double rms_deltas() const
    {
      return std::sqrt(deltas_.dot(deltas_)/9);
    }

    //! The gradient of R = w(sum(delta))^2 with respect to U_cart
    scitbx::sym_mat3<double>
    gradients() const
    {
      scitbx::sym_mat3<double> gradients;
      for (int i=0;i<3;i++) {
        gradients[i] = weight * 2 * deltas_[i];
      }
      for (int i=3;i<6;i++) {
        gradients[i] = weight * 4 * deltas_[i];
      }
      return gradients;
    }

    void
    linearise(
      cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
      cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
      unsigned const& i_seq) const
    {
      cctbx::xray::parameter_indices const &ids
        = parameter_map[i_seq];
      // One restraint per parameter == six rows in the restraint matrix
      for (std::size_t i=0;i<6;i++) {
        std::size_t row_i = linearised_eqns.next_row();
        if (i < 3) {
          for (std::size_t j=0;j<3;j++) {
            if (j==i) {
              linearised_eqns.design_matrix(row_i, ids.u_aniso+j) = 2./3;
            }
            else {
              linearised_eqns.design_matrix(row_i, ids.u_aniso+j) = -1./3;
            }
          }
        }
        else {
          linearised_eqns.design_matrix(row_i, ids.u_aniso+i) = 2.;
        }
        linearised_eqns.weights[row_i] = weight;
        linearised_eqns.deltas[row_i] = deltas_[i];
      }
    }

    //! Support for isotropic_adp_residual_sum.
    /*! Not available in Python.
     */
    void
    add_gradients(
      af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart,
      unsigned const& i_seq) const
    {
      gradients_aniso_cart[i_seq] += gradients();
    }

    //! Cartesian anisotropic displacement parameters.
    scitbx::sym_mat3<double> u_cart;
    double weight;

  protected:
    scitbx::sym_mat3<double> deltas_;

    void init_deltas()
    {
      double const u_iso =
        adptbx::u_cart_as_u_iso(u_cart);
      for (int i=0;i<3;i++) {
        deltas_[i] = u_cart[i] - u_iso;
      }
      for (int i=3;i<6;i++) {
        deltas_[i] = u_cart[i];
      }
    }

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
    af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart)
  {
    CCTBX_ASSERT(   gradients_aniso_cart.size() == 0
                 || gradients_aniso_cart.size() == u_cart.size());
    double result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      isotropic_adp_proxy const& proxy = proxies[i];
      isotropic_adp restraint(u_cart, proxy);
      result += restraint.residual();
      if (gradients_aniso_cart.size() != 0) {
        restraint.add_gradients(gradients_aniso_cart, proxy.i_seq);
      }
    }
    return result;
  }

  /*! \brief Fast computation of isotropic_adp::residual() given an array
      of adp_similarity proxies.
   */
  af::shared<double>
  isotropic_adp_residuals(
    af::const_ref<scitbx::sym_mat3<double> > const& u_cart,
    af::const_ref<isotropic_adp_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      isotropic_adp_proxy const& proxy = proxies[i];
      isotropic_adp restraint(u_cart, proxy);
      result.push_back(restraint.residual());
    }
    return result;
  }

  /*! \brief Fast computation of isotropic_adp::rms_deltas() given an array
      of isotropic_adp proxies.
   */
  af::shared<double>
  isotropic_adp_deltas_rms(
    af::const_ref<scitbx::sym_mat3<double> > const& u_cart,
    af::const_ref<isotropic_adp_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      result.push_back(isotropic_adp(u_cart, proxies[i]).rms_deltas());
    }
    return result;
  }

}} // cctbx::adp_restraints

#endif
