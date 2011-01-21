#ifndef CCTBX_ADP_RESTRAINTS_ADP_SIMILARITY_H
#define CCTBX_ADP_RESTRAINTS_ADP_SIMILARITY_H

#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <cctbx/adptbx.h>
#include <cctbx/restraints.h>
#include <scitbx/matrix/matrix_vector_operations.h>

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
      adp_similarity_proxy const& proxy)
    :
      weight(proxy.weight)
    {
      CCTBX_ASSERT(u_cart_.size() == u_iso_.size());
      CCTBX_ASSERT(u_cart_.size() == use_u_aniso_.size());
      for (int i=0;i<2;i++) {
        std::size_t i_seq = proxy.i_seqs[i];
        CCTBX_ASSERT(i_seq < u_cart_.size());
        u_cart[i] = u_cart_[i_seq];
        u_iso[i] = u_iso_[i_seq];
        use_u_aniso[i] = use_u_aniso_[i_seq];
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

    scitbx::sym_mat3<double>
    gradient_0() const
    {
      scitbx::sym_mat3<double> result;
      for(int i=0;i<3;i++) {
        result[i] = 2 * deltas_[i];
      }
      for(int i=3;i<6;i++) {
        result[i] = 4 * deltas_[i];
      }
      return result;
    }

    //! This returns gradients_u_cart and gradients_u_equiv combined
    af::tiny<scitbx::sym_mat3<double>, 2>
    gradients() const
    {
      af::tiny<scitbx::sym_mat3<double>, 2> result;
      result[0] = gradient_0();
      result[1] = -result[0];
      return result;
    }


    void
    linearise(
      uctbx::unit_cell const &unit_cell,
      cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
      cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
      af::tiny<unsigned, 2> const& i_seqs) const
    {
      static const double grads_u_cart[6][6] = {
        { 1, 0, 0, 0, 0, 0},
        { 0, 1, 0, 0, 0, 0},
        { 0, 0, 1, 0, 0, 0},
        { 0, 0, 0, 1, 0, 0},
        { 0, 0, 0, 0, 1, 0},
        { 0, 0, 0, 0, 0, 1},
      };
      if (!use_u_aniso[0] && !use_u_aniso[1]) {
        // Only one restraint, i.e. one row added to restraint matrix
        std::size_t row_i = linearised_eqns.next_row();
        for (std::size_t j=0;j<2;j++) {
          double grad = 1.0;
          if (j == 1) grad *= -1;
          cctbx::xray::parameter_indices const &ids_j
            = parameter_map[i_seqs[j]];
          if (ids_j.u_iso == -1) continue;
          linearised_eqns.design_matrix(row_i, ids_j.u_iso) = grad;
        }
        linearised_eqns.weights[row_i] = weight;
        linearised_eqns.deltas[row_i] = deltas_[0];
      }
      else {
        // One restraint per parameter == six rows in the restraint matrix
        af::const_ref<double, af::mat_grid> const &f
          = unit_cell.u_star_to_u_cart_linear_map();
        for (std::size_t i=0;i<6;i++) {
          scitbx::sym_mat3<double> grad_u_star;
          scitbx::matrix::matrix_transposed_vector(
            6, 6, f.begin(),
            scitbx::sym_mat3<double>(grads_u_cart[i]).begin(),
            grad_u_star.begin());
          std::size_t row_i = linearised_eqns.next_row();
          for (std::size_t j=0;j<2;j++) {
            if (j == 1) {
              grad_u_star = -grad_u_star;
            }
            cctbx::xray::parameter_indices const &ids_j
              = parameter_map[i_seqs[j]];
            if (use_u_aniso[j] && ids_j.u_aniso != -1) {
              for (std::size_t k=0;k<6;k++) {
                double grad_k = grad_u_star[k];
                if (k>2) grad_k *= 2; // symmetric matrix, off-diagonals count double
                linearised_eqns.design_matrix(row_i, ids_j.u_aniso+k) = grad_k;
              }
            }
            else if (i < 3 && !use_u_aniso[j] && ids_j.u_iso != -1) {
              for (std::size_t k=0;k<3;k++) {
                linearised_eqns.design_matrix(row_i, ids_j.u_iso)
                  = grad_u_star[k];
              }
            }
          }
          linearised_eqns.weights[row_i] = weight;
          linearised_eqns.deltas[row_i] = deltas_[i];
        }
      }
    }

    //! Support for adp_similarity_residual_sum.
    /*! Not available in Python.
     */
    void
    add_gradients(
      af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart,
      af::ref<double> const& gradients_iso,
      af::tiny<unsigned, 2> const& i_seqs) const
    {
      //! () - ()
      if (use_u_aniso[0] && use_u_aniso[1]) {
        scitbx::sym_mat3<double> g0 = gradient_0();
        gradients_aniso_cart[i_seqs[0]] += g0;
        gradients_aniso_cart[i_seqs[1]] += -g0;
      }
      //! () - o
      else if (use_u_aniso[0] && !use_u_aniso[1]) {
        scitbx::sym_mat3<double> g0 = gradient_0();
        gradients_aniso_cart[i_seqs[0]] += g0;
        gradients_iso[i_seqs[1]] += -g0.trace();
      }
      //! o - ()
      else if (!use_u_aniso[0] && use_u_aniso[1]) {
        scitbx::sym_mat3<double> g0 = gradient_0();
        gradients_iso[i_seqs[0]] += g0.trace();
        gradients_aniso_cart[i_seqs[1]] += -g0;
      }
      //! o - o
      else if (!use_u_aniso[0] && !use_u_aniso[1]) {
        double g_iso = 2 * deltas_[0];
        gradients_iso[i_seqs[0]] += g_iso;
        gradients_iso[i_seqs[1]] += -g_iso;
      }
    }

    //! Cartesian anisotropic displacement parameters.
    af::tiny<scitbx::sym_mat3<double>, 2> u_cart;
    //! Isotropic displacement parameters.
    af::tiny<double, 2> u_iso;
    //! Use U aniso.
    af::tiny<bool, 2> use_u_aniso;

    double weight;

  protected:

    scitbx::sym_mat3<double> deltas_;

    void init_deltas()
    {
      deltas_ = scitbx::sym_mat3<double> (0,0,0,0,0,0);
      //! () - ()
      if (use_u_aniso[0] && use_u_aniso[1]) {
        scitbx::sym_mat3<double> const& u_1 = u_cart[0];
        scitbx::sym_mat3<double> const& u_2 = u_cart[1];
        for (int i=0;i<6;i++) {
          deltas_[i] = u_1[i] - u_2[i];
        }
      }
      //! () - o
      else if (use_u_aniso[0] && !use_u_aniso[1]) {
        scitbx::sym_mat3<double> const& u_1 = u_cart[0];
        double const& u_2 = u_iso[1];
        scitbx::sym_mat3<double> const& u_2_cart = adptbx::u_iso_as_u_cart(u_2);
        for (int i=0;i<6;i++) {
          deltas_[i] = u_1[i] - u_2_cart[i];
        }
      }
      //! o - ()
      else if (!use_u_aniso[0] && use_u_aniso[1]) {
        double const& u_1 = u_iso[0];
        scitbx::sym_mat3<double> const& u_2 = u_cart[1];
        scitbx::sym_mat3<double> const& u_1_cart = adptbx::u_iso_as_u_cart(u_1);
        for (int i=0;i<6;i++) {
          deltas_[i] = u_1_cart[i] - u_2[i];
        }
      }
      //! o - o
      else if (!use_u_aniso[0] && !use_u_aniso[1]) {
        double const& u_1 = u_iso[0];
        double const& u_2 = u_iso[1];
        deltas_[0] = u_1 - u_2;
      }
    }

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
    af::ref<double> const& gradients_iso)
  {
    CCTBX_ASSERT(   gradients_aniso_cart.size() == 0
                 || gradients_aniso_cart.size() == u_cart.size());
    CCTBX_ASSERT(gradients_aniso_cart.size() == gradients_iso.size());
    double result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      adp_similarity_proxy const& proxy = proxies[i];
      adp_similarity restraint(u_cart, u_iso, use_u_aniso, proxy);
      result += restraint.residual();
      if (gradients_aniso_cart.size() != 0) {
        restraint.add_gradients(gradients_aniso_cart, gradients_iso, proxy.i_seqs);
      }
    }
    return result;
  }

  /*! \brief Fast computation of adp_similarity::residual() given an array
      of adp_similarity proxies.
   */
  af::shared<double>
  adp_similarity_residuals(
    af::const_ref<scitbx::sym_mat3<double> > const& u_cart,
    af::const_ref<double> const& u_iso,
    af::const_ref<bool> const& use_u_aniso,
    af::const_ref<adp_similarity_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      adp_similarity_proxy const& proxy = proxies[i];
      adp_similarity restraint(u_cart, u_iso, use_u_aniso, proxy);
      result.push_back(restraint.residual());
    }
    return result;
  }

  /*! \brief Fast computation of adp_similarity::rms_deltas() given an array
      of adp_similarity proxies.
   */
  af::shared<double>
  adp_similarity_deltas_rms(
    af::const_ref<scitbx::sym_mat3<double> > const& u_cart,
    af::const_ref<double> const& u_iso,
    af::const_ref<bool> const& use_u_aniso,
    af::const_ref<adp_similarity_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      result.push_back(adp_similarity(u_cart, u_iso, use_u_aniso, proxies[i]).rms_deltas());
    }
    return result;
  }

}} // cctbx::adp_restraints

#endif
