#ifndef CCTBX_ADP_RESTRAINTS_H
#define CCTBX_ADP_RESTRAINTS_H

#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <cctbx/adptbx.h>
#include <cctbx/restraints.h>
#include <scitbx/matrix/matrix_vector_operations.h>

namespace cctbx { namespace adp_restraints {

  template <typename FloatType>
  struct adp_restraint_params {
    af::shared<scitbx::vec3<FloatType> > sites_cart;
    af::shared<scitbx::sym_mat3<FloatType> > u_cart;
    af::shared<FloatType> u_iso;
    af::shared<bool> use_u_aniso;

    adp_restraint_params(
      af::shared<scitbx::vec3<FloatType> > const &sites_cart_,
      af::shared<scitbx::sym_mat3<FloatType> > const &u_cart_,
      af::shared<FloatType> const &u_iso_,
      af::shared<bool> const &use_u_aniso_)
      : sites_cart(sites_cart_),
        u_cart(u_cart_),
        u_iso(u_iso_),
        use_u_aniso(use_u_aniso_)
    {}

    adp_restraint_params(
      af::shared<scitbx::sym_mat3<FloatType> > const &u_cart_,
      af::shared<FloatType> const &u_iso_,
      af::shared<bool> const &use_u_aniso_)
      : u_cart(u_cart_),
        u_iso(u_iso_),
        use_u_aniso(use_u_aniso_)
    {}

    adp_restraint_params(
      af::shared<scitbx::vec3<FloatType> > const &sites_cart_,
      af::shared<scitbx::sym_mat3<FloatType> > const &u_cart_)
      : sites_cart(sites_cart_),
        u_cart(u_cart_),
        use_u_aniso(u_cart_.size())
    {
      for (int i=0; i < use_u_aniso.size(); i++)
        use_u_aniso[i] = true;
    }

    adp_restraint_params(
      af::shared<scitbx::sym_mat3<FloatType> > const &u_cart_)
      : u_cart(u_cart_),
        use_u_aniso(u_cart_.size())
    {
      for (int i=0; i < use_u_aniso.size(); i++)
        use_u_aniso[i] = true;
    }

    adp_restraint_params(
      af::shared<FloatType> const &u_iso_)
      : u_iso(u_iso_),
        use_u_aniso(u_iso_.size())
    {
      for (int i=0; i < use_u_aniso.size(); i++)
        use_u_aniso[i] = false;
    }

  };

  template <typename GradientSource>
  void linearise_1(uctbx::unit_cell const &unit_cell,
    cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
    cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
    unsigned i_seq,
    bool use_u_aniso,
    double weight,
    scitbx::sym_mat3<double> const &deltas)
  {
    cctbx::xray::parameter_indices const &ids = parameter_map[i_seq];
    if (use_u_aniso) {
      // One restraint per parameter == six rows in the restraint matrix
      CCTBX_ASSERT(ids.u_aniso != -1);
      for (int i=0; i<6; i++) {
        std::size_t row_i = linearised_eqns.next_row();
        scitbx::sym_mat3<double> grad_u_star;
        scitbx::matrix::matrix_transposed_vector(
          6, 6, unit_cell.u_star_to_u_cart_linear_map().begin(),
          scitbx::sym_mat3<double>(GradientSource::cart_grad_row(i)).begin(),
          grad_u_star.begin());
        for (int j=0; j<6; j++) {
          // symmetric matrix, off-diagonals count double
          linearised_eqns.design_matrix(row_i, ids.u_aniso+j) =
            (j > 2 ? 2*grad_u_star[j] : grad_u_star[j]);
        }
        linearised_eqns.weights[row_i] = weight;
        linearised_eqns.deltas[row_i] = deltas[i];
      }
    }
    else {
      CCTBX_ASSERT(ids.u_iso != -1);
      std::size_t row_i = linearised_eqns.next_row();
      linearised_eqns.design_matrix(row_i, ids.u_iso) =
        GradientSource::grad_u_iso(0);
      linearised_eqns.weights[row_i] = weight;
      linearised_eqns.deltas[row_i] = deltas[0];
    }
  }

  template <typename GradientSource>
  void
    linearise_2(
    GradientSource const &grad_src,
    uctbx::unit_cell const &unit_cell,
    cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
    cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
    af::tiny<unsigned, 2> const &i_seqs,
    af::tiny<bool, 2> const &use_u_aniso,
    double weight,
    scitbx::sym_mat3<double> const &deltas)
  {
    if (!use_u_aniso[0] && !use_u_aniso[1]) {
      // Only one restraint, i.e. one row added to restraint matrix
      std::size_t row_i = linearised_eqns.next_row();
      for (int j=0; j<2; j++) {
        cctbx::xray::parameter_indices const &ids_j =
          parameter_map[i_seqs[j]];
        if (ids_j.u_iso == -1) continue;
        linearised_eqns.design_matrix(row_i, ids_j.u_iso) =
          (j == 0 ? 1 : -1) * grad_src.grad_u_iso(j);
      }
      linearised_eqns.weights[row_i] = weight;
      linearised_eqns.deltas[row_i] = deltas[0];
    }
    else if(use_u_aniso[0] && use_u_aniso[1]) {
      scitbx::sym_mat3<double> grad_u_star[2];
      for (int i=0; i<6; i++) {
        std::size_t row_i = linearised_eqns.next_row();
        scitbx::matrix::matrix_transposed_vector(
          6, 6, unit_cell.u_star_to_u_cart_linear_map().begin(),
          scitbx::sym_mat3<double>(grad_src.cart_grad_row(0, i)).begin(),
          grad_u_star[0].begin());
        scitbx::matrix::matrix_transposed_vector(
          6, 6, unit_cell.u_star_to_u_cart_linear_map().begin(),
          scitbx::sym_mat3<double>(grad_src.cart_grad_row(1, i)).begin(),
          grad_u_star[1].begin());
        grad_u_star[1] = -grad_u_star[1];
        for (int j=0; j<2; j++) {
          cctbx::xray::parameter_indices const &ids_j
            = parameter_map[i_seqs[j]];
          if (use_u_aniso[j] && ids_j.u_aniso != -1) {
            for (int k=0; k<6; k++) {
              linearised_eqns.design_matrix(row_i, ids_j.u_aniso+k) =
                grad_u_star[j][k] * (k > 2 ? 2 : 1);
            }
          }
        }
        linearised_eqns.weights[row_i] = weight;
        linearised_eqns.deltas[row_i] = deltas[i];
      }
    }
    else {
      scitbx::sym_mat3<double> grad_u_star;
      double grad_u_iso;
      for (int i=0; i<6; i++) {
        std::size_t row_i = linearised_eqns.next_row();
        if (use_u_aniso[0]) {
          scitbx::matrix::matrix_transposed_vector(
            6, 6, unit_cell.u_star_to_u_cart_linear_map().begin(),
            scitbx::sym_mat3<double>(grad_src.cart_grad_row(0, i)).begin(),
            grad_u_star.begin());
          grad_u_iso = -grad_src.grad_u_iso(1);
        }
        else {
          scitbx::matrix::matrix_transposed_vector(
            6, 6, unit_cell.u_star_to_u_cart_linear_map().begin(),
            scitbx::sym_mat3<double>(grad_src.cart_grad_row(1, i)).begin(),
            grad_u_star.begin());
          grad_u_star = -grad_u_star;
          grad_u_iso = grad_src.grad_u_iso(0);
        }
        for (int j=0; j<2; j++) {
          cctbx::xray::parameter_indices const &ids_j
            = parameter_map[i_seqs[j]];
          if (use_u_aniso[j] && ids_j.u_aniso != -1) {
            for (int k=0; k<6; k++) {
              linearised_eqns.design_matrix(row_i, ids_j.u_aniso+k) =
                grad_u_star[k] * (k > 2 ? 2 : 1);
            }
          }
          else if (i < 3 && !use_u_aniso[j] && ids_j.u_iso != -1) {
            linearised_eqns.design_matrix(row_i, ids_j.u_iso) =
              grad_u_iso;
          }
        }
        linearised_eqns.weights[row_i] = weight;
        linearised_eqns.deltas[row_i] = deltas[i];
      }
    }
  }

  /* as above, but optimised for the static gradients */
  template <typename GradientSource>
  void
    linearise_2(
    uctbx::unit_cell const &unit_cell,
    cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
    cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
    af::tiny<unsigned, 2> const &i_seqs,
    af::tiny<bool, 2> const &use_u_aniso,
    double weight,
    scitbx::sym_mat3<double> const &deltas)
  {
    if (!use_u_aniso[0] && !use_u_aniso[1]) {
      // Only one restraint, i.e. one row added to restraint matrix
      std::size_t row_i = linearised_eqns.next_row();
      for (int j=0; j<2; j++) {
        cctbx::xray::parameter_indices const &ids_j =
          parameter_map[i_seqs[j]];
        if (ids_j.u_iso == -1) continue;
        linearised_eqns.design_matrix(row_i, ids_j.u_iso) =
          (j == 0 ? 1 : -1) * GradientSource::grad_u_iso(j);
      }
      linearised_eqns.weights[row_i] = weight;
      linearised_eqns.deltas[row_i] = deltas[0];
    }
    else {
      scitbx::sym_mat3<double> grad_u_star;
      for (int i=0; i<6; i++) {
        scitbx::matrix::matrix_transposed_vector(
          6, 6, unit_cell.u_star_to_u_cart_linear_map().begin(),
          scitbx::sym_mat3<double>(GradientSource::cart_grad_row(i)).begin(),
          grad_u_star.begin());
        std::size_t row_i = linearised_eqns.next_row();
        for (int j=0; j<2; j++) {
          if (j == 1)
            grad_u_star = -grad_u_star;
          cctbx::xray::parameter_indices const &ids_j
            = parameter_map[i_seqs[j]];
          if (use_u_aniso[j] && ids_j.u_aniso != -1) {
            for (int k=0; k<6; k++) {
              linearised_eqns.design_matrix(row_i, ids_j.u_aniso+k) =
                grad_u_star[k] * (k > 2 ? 2 : 1);
            }
          }
          else if (i < 3 && !use_u_aniso[j] && ids_j.u_iso != -1) {
            linearised_eqns.design_matrix(row_i, ids_j.u_iso) =
              (j == 0 ? 1 : -1) * GradientSource::grad_u_iso(j);
          }
        }
        linearised_eqns.weights[row_i] = weight;
        linearised_eqns.deltas[row_i] = deltas[i];
      }
    }
  }

  template <int n_adp> struct adp_restraint_proxy {
    adp_restraint_proxy() {}

    //! Constructor.
    adp_restraint_proxy(af::tiny<unsigned, n_adp> const &i_seqs_,
      double weight_)
    : i_seqs(i_seqs_),
      weight(weight_)
    {}

    //! Indices into array of sites.
    af::tiny<unsigned, n_adp> i_seqs;
    //! weight
    double weight;
  };


  template <int n_adp> class adp_restraint_base {
  public:

    adp_restraint_base() : weight(0)  {}

    adp_restraint_base(double weight_) : weight(weight_) {}

    adp_restraint_base(af::tiny<bool, n_adp> const &use_u_aniso_,
      double weight_)
    : use_u_aniso(use_u_aniso_),
      weight(weight_)
    {}

    adp_restraint_base(
      adp_restraint_params<double> const &params,
      adp_restraint_proxy<n_adp> const &proxy)
    : weight(proxy.weight)
    {
      for (int i=0; i<n_adp; i++) {
        std::size_t i_seq = proxy.i_seqs[i];
        CCTBX_ASSERT(i_seq < params.use_u_aniso.size());
        use_u_aniso[i] = params.use_u_aniso[i_seq];
      }
    }

    scitbx::sym_mat3<double> deltas() const { return deltas_; }

    //! weight * [[sum_{ii} (deltas)**2] + [2 * sum_{i<j} (deltas)**2]].
    /* This is the square of the Frobenius norm of the matrix of deltas, or
       alternatively the inner product of the matrix of deltas with itself.
       This is used since the residual is then rotationally invariant.
     */
    double residual() const { return weight * deltas_.dot(deltas_); }

    //! sqrt(mean_sq(deltas))
    //! The off-diagonal elements are included twice.
    double rms_deltas() const { return std::sqrt(deltas_.dot(deltas_)/9); }

    //! The gradient of R = w(sum(delta))^2 with respect to U_cart
    scitbx::sym_mat3<double> gradients() const {
      scitbx::sym_mat3<double> gradients;
      for (int i=0; i<6; i++)
        gradients[i] = weight*deltas_[i]*(i < 3 ? 2 : 4);
      return gradients;
    }

    //! Use U aniso.
    af::tiny<bool, n_adp> use_u_aniso;
    // restraint weight
    double weight;
  protected:
    scitbx::sym_mat3<double> deltas_;
  };

  /*! Fast computation of sum of fixed_u_eq_adp::residual() and gradients
      given an array of isotropic_adp proxies.
   */
  /*! The fixed_u_eq_adp::gradients() are added to the gradient_array if
      gradient_array.size() == sites_cart.size().
      gradient_array must be initialized before this function
      is called.
      No gradient calculations are performed if gradient_array.size() == 0.
   */
  template <typename ProxyType, typename RestraintType>
  double adp_restraint_residual_sum(
    adp_restraint_params<double> const &params,
    af::const_ref<ProxyType> const& proxies,
    af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart,
    af::ref<double> const& gradients_iso)
  {
    CCTBX_ASSERT(gradients_aniso_cart.size() == 0 ||
      gradients_aniso_cart.size() == params.u_cart.size());
    CCTBX_ASSERT(gradients_aniso_cart.size() == gradients_iso.size());
    double result = 0;
    for(std::size_t i=0; i<proxies.size(); i++) {
      RestraintType restraint(params, proxies[i]);
      result += restraint.residual();
      if (gradients_aniso_cart.size() != 0) {
        restraint.add_gradients(
          gradients_aniso_cart, gradients_iso, proxies[i].i_seqs);
      }
    }
    return result;
  }

  /* similar to the function above - specialised for anisotropic gradients only
  */
  template <typename ProxyType, typename RestraintType>
  double adp_restraint_residual_sum_aniso(
    adp_restraint_params<double> const &params,
    af::const_ref<ProxyType> const& proxies,
    af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart)
  {
    CCTBX_ASSERT(gradients_aniso_cart.size() == 0 ||
      gradients_aniso_cart.size() == params.u_cart.size());
    double result = 0;
    for(std::size_t i=0; i<proxies.size(); i++) {
      RestraintType restraint(params, proxies[i]);
      result += restraint.residual();
      if (gradients_aniso_cart.size() != 0) {
        restraint.add_gradients(gradients_aniso_cart, proxies[i].i_seqs);
      }
    }
    return result;
  }

  /*! \brief Fast computation of isotropic_adp::residual() given an array
      of the proxies.
   */
  template <typename ProxyType, typename RestraintType>
  af::shared<double> adp_restraint_residuals(
    adp_restraint_params<double> const &params,
    af::const_ref<ProxyType> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0; i<proxies.size(); i++) {
      result.push_back(
        RestraintType(params, proxies[i]).residual());
    }
    return result;
  }

  /*! \brief Fast computation of fixed_u_eq_adp::rms_deltas() given an array
      of the proxies.
   */
  template <typename ProxyType, typename RestraintType>
  af::shared<double> adp_restraint_deltas_rms(
    adp_restraint_params<double> const &params,
    af::const_ref<ProxyType> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0; i<proxies.size(); i++) {
      result.push_back(
        RestraintType(params, proxies[i]).rms_deltas());
    }
    return result;
  }

}} // cctbx::adp_restraints
#endif
