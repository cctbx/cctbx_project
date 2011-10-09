#ifndef CCTBX_ADP_RESTRAINTS_ADP_SIMILARITY_H
#define CCTBX_ADP_RESTRAINTS_ADP_SIMILARITY_H

#include <cctbx/adp_restraints/adp_restraints.h>

namespace cctbx { namespace adp_restraints {

  using scitbx::vec3;
  using scitbx::mat3;

  struct adp_similarity_proxy : public adp_restraint_proxy<2> {
    adp_similarity_proxy() {}
    adp_similarity_proxy(
      af::tiny<unsigned, 2> const& i_seqs,
      double weight)
    : adp_restraint_proxy<2>(i_seqs, weight)
    {}
  };

  class adp_similarity : public adp_restraint_base_6<2> {
  public:
    //! Constructor.
    adp_similarity(
      af::tiny<scitbx::sym_mat3<double>, 2> const& u_cart,
      double weight)
    : adp_restraint_base_6<2>(af::tiny<bool, 2>(true, true), weight)
    {
      init_deltas(u_cart[0], u_cart[1]);
    }

    adp_similarity(
      af::tiny<double, 2> const& u_iso,
      double weight)
    : adp_restraint_base_6<2>(af::tiny<bool, 2>(false, false), weight)
    {
      init_deltas(u_iso[0], u_iso[1]);
    }

    adp_similarity(
      scitbx::sym_mat3<double> const& u_cart,
      double u_iso,
      double weight)
    : adp_restraint_base_6<2>(af::tiny<bool, 2>(true, false), weight)
    {
      init_deltas(u_cart, u_iso);
    }

    adp_similarity(
      double u_iso,
      scitbx::sym_mat3<double> const& u_cart,
      double weight)
    : adp_restraint_base_6<2>(af::tiny<bool, 2>(false, true), weight)
    {
      init_deltas(u_iso, u_cart);
    }

    //! Constructor.
    adp_similarity(
      adp_restraint_params<double> const &params,
      adp_similarity_proxy const& proxy)
    : adp_restraint_base_6<2>(params, proxy)
    {
      if (use_u_aniso[0] && use_u_aniso[1]) {
        CCTBX_ASSERT(proxy.i_seqs[0] < params.u_cart.size());
        CCTBX_ASSERT(proxy.i_seqs[1] < params.u_cart.size());
        init_deltas(params.u_cart[proxy.i_seqs[0]], params.u_cart[proxy.i_seqs[1]]);
      }
      else if (use_u_aniso[0] && !use_u_aniso[1]) {
        CCTBX_ASSERT(proxy.i_seqs[0] < params.u_cart.size());
        CCTBX_ASSERT(proxy.i_seqs[1] < params.u_iso.size());
        init_deltas(params.u_cart[proxy.i_seqs[0]], params.u_iso[proxy.i_seqs[1]]);
      }
      else if (!use_u_aniso[0] && use_u_aniso[1]) {
        CCTBX_ASSERT(proxy.i_seqs[0] < params.u_iso.size());
        CCTBX_ASSERT(proxy.i_seqs[1] < params.u_cart.size());
        init_deltas(params.u_iso[proxy.i_seqs[0]], params.u_cart[proxy.i_seqs[1]]);
      }
      else {
        CCTBX_ASSERT(proxy.i_seqs[0] < params.u_iso.size());
        CCTBX_ASSERT(proxy.i_seqs[1] < params.u_iso.size());
        init_deltas(params.u_iso[proxy.i_seqs[0]], params.u_iso[proxy.i_seqs[1]]);
      }
    }

    void linearise(
      uctbx::unit_cell const &unit_cell,
      cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
      cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
      af::tiny<unsigned, 2> const& i_seqs) const
    {
      linearise_2<adp_similarity>(
        unit_cell, linearised_eqns, parameter_map, i_seqs, use_u_aniso, weight,
        deltas_.begin());
    }

    static double grad_u_iso(int) {  return 1;  }

    static const double* cart_grad_row(int i) {
      static const double grads_u_cart[6][6] = {
        { 1, 0, 0, 0, 0, 0},
        { 0, 1, 0, 0, 0, 0},
        { 0, 0, 1, 0, 0, 0},
        { 0, 0, 0, 1, 0, 0},
        { 0, 0, 0, 0, 1, 0},
        { 0, 0, 0, 0, 0, 1},
      };
      return &grads_u_cart[i][0];
    }

  protected:

    void init_deltas(scitbx::sym_mat3<double> const &u_cart1,
      scitbx::sym_mat3<double> const &u_cart2)
    {
      deltas_ = u_cart1 - u_cart2;
    }

    void init_deltas(double u_iso1, double u_iso2) {
      deltas_[0] = u_iso1 - u_iso2;
      for (int i=1; i<6; i++) deltas_[i] = 0;
    }

    void init_deltas(scitbx::sym_mat3<double> const &u_cart, double u_iso) {
      for (int i=0; i<6; i++)
        deltas_[i] = u_cart[i] - (i < 3 ? u_iso : 0);
    }

    void init_deltas(double u_iso, scitbx::sym_mat3<double> const &u_cart) {
      for (int i=0; i<6; i++)
        deltas_[i] = (i < 3 ? u_iso : 0 ) - u_cart[i];
    }

  };


struct adp_u_eq_similarity_proxy : public adp_restraint_proxy_n {
    adp_u_eq_similarity_proxy() {}
    adp_u_eq_similarity_proxy(
      af::shared<unsigned> const& i_seqs,
      double weight)
    : adp_restraint_proxy_n(i_seqs, weight)
    {}
  };

class adp_u_eq_similarity : public adp_restraint_base_n {
  public:
    //! Constructor.
    adp_u_eq_similarity(
      adp_restraint_params<double> const &params,
      adp_u_eq_similarity_proxy const& proxy)
    : adp_restraint_base_n(params, proxy),
      mean_u_eq(0)
    {
      for (int i=0; i < proxy.i_seqs.size(); i++) {
        if (use_u_aniso[i]) {
          CCTBX_ASSERT(proxy.i_seqs[i] < params.u_cart.size());
          deltas_[i] = params.u_cart[proxy.i_seqs[i]].trace()/3;
          mean_u_eq += deltas_[i];
        }
        else {
          CCTBX_ASSERT(proxy.i_seqs[i] < params.u_iso.size());
          deltas_[i] = params.u_iso[proxy.i_seqs[i]];
          mean_u_eq += deltas_[i];
        }
      }
      mean_u_eq /= proxy.i_seqs.size();
      for (int i=0; i < proxy.i_seqs.size(); i++)
        deltas_[i] -= mean_u_eq;
    }

    void linearise(uctbx::unit_cell const &unit_cell,
      cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
      cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
      af::shared<unsigned> const &i_seqs) const
    {
      CCTBX_ASSERT(use_u_aniso.size()==i_seqs.size());
      double k_ij = -1./(3*deltas_.size()), k_ii = 1./3 + k_ij;
      scitbx::sym_mat3<double> u_star_grad_ii, u_star_grad_ij;
      scitbx::matrix::matrix_transposed_vector(
        6, 6, unit_cell.u_star_to_u_cart_linear_map().begin(),
        (scitbx::sym_mat3<double>(k_ii)).begin(),
        u_star_grad_ii.begin());
      scitbx::matrix::matrix_transposed_vector(
        6, 6, unit_cell.u_star_to_u_cart_linear_map().begin(),
        (scitbx::sym_mat3<double>(k_ij)).begin(),
        u_star_grad_ij.begin());
      double u_iso_ij = -1./deltas_.size(), u_iso_ii = 1 + u_iso_ij;
      for (int i=0; i < i_seqs.size(); i++) {
        std::size_t row_i = linearised_eqns.next_row();
        for (int j=0; j < i_seqs.size(); j++) {
          cctbx::xray::parameter_indices const &jds = parameter_map[i_seqs[j]];
          if (use_u_aniso[j]) {
            CCTBX_ASSERT(jds.u_aniso != -1);
            scitbx::sym_mat3<double> &u_cart_grad =
              (i == j ? u_star_grad_ii : u_star_grad_ij);
            for (int k=0; k < 6; k++) {
              linearised_eqns.design_matrix(row_i, jds.u_aniso+k) =
                (k > 2 ? 2*u_cart_grad[k] : u_cart_grad[k]);
            }
          }
          else {
            CCTBX_ASSERT(jds.u_iso != -1);
            linearised_eqns.design_matrix(row_i, jds.u_iso) =
              (i == j ? u_iso_ii : u_iso_ij);
          }
        }
        linearised_eqns.weights[row_i] = weight;
        linearised_eqns.deltas[row_i] = deltas_[i];
      }
    }

    double mean_u_eq;
  };

  struct adp_volume_similarity_proxy : public adp_restraint_proxy_n {
    adp_volume_similarity_proxy() {}
    adp_volume_similarity_proxy(
      af::shared<unsigned> const& i_seqs,
      double weight)
    : adp_restraint_proxy_n(i_seqs, weight)
    {}
  };

/* in this restraint the gradients are estimated considering that eigen
values and eigen vectors are independent */
  class adp_volume_similarity : public adp_restraint_base_n {
  public:

    adp_volume_similarity(
      adp_restraint_params<double> const &params,
      adp_volume_similarity_proxy const& proxy)
    : adp_restraint_base_n(params, proxy),
      mean_u_volume(0),
      grad_indices(proxy.i_seqs.size())
    {
      std::size_t u_iso_idx = 0, u_star_idx = 0;
      for (int i=0; i < proxy.i_seqs.size(); i++) {
        if (use_u_aniso[i]) {
          CCTBX_ASSERT(proxy.i_seqs[i] < params.u_cart.size());
          adptbx::eigensystem<double> es(params.u_cart[proxy.i_seqs[i]]);
          const vec3<double> &v = es.values();
          deltas_[i] = std::sqrt(v[0]*v[1]*v[2]);
          mean_u_volume += deltas_[i];
          u_cart_grads.push_back(calc_grad(es));
          grad_indices[i] = u_star_idx++;
        }
        else {
          CCTBX_ASSERT(proxy.i_seqs[i] < params.u_iso.size());
          deltas_[i] = std::pow(params.u_iso[proxy.i_seqs[i]], 3./2);
          mean_u_volume += deltas_[i];
          iso_grads.push_back(
            scitbx::constants::two_pi*std::sqrt(
              params.u_iso[proxy.i_seqs[i]]));
          grad_indices[i] = u_iso_idx++;
        }
      }
      mean_u_volume /= proxy.i_seqs.size();
      for (int i=0; i < proxy.i_seqs.size(); i++)
        deltas_[i] = r3diff_to_vol(deltas_[i]-mean_u_volume);
      mean_u_volume = r3diff_to_vol(mean_u_volume);
    }

    void linearise(uctbx::unit_cell const &unit_cell,
      cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
      cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
      af::shared<unsigned> const &i_seqs) const
    {
      CCTBX_ASSERT(use_u_aniso.size()==i_seqs.size());
      double k_ij = -1./deltas_.size(), k_ii = 1+k_ij;
      std::size_t first_row_i = linearised_eqns.next_row();
      for (int i=0; i < i_seqs.size(); i++) {
        if (i != 0)
          linearised_eqns.next_row();
        cctbx::xray::parameter_indices const &ids = parameter_map[i_seqs[i]];
        if (use_u_aniso[i]) {
          CCTBX_ASSERT(ids.u_aniso != -1);
          scitbx::sym_mat3<double> u_star_grad_ii, u_star_grad_ij;
          scitbx::matrix::matrix_transposed_vector(
            6, 6, unit_cell.u_star_to_u_cart_linear_map().begin(),
            (u_cart_grads[grad_indices[i]]*k_ii).begin(),
            u_star_grad_ii.begin());
          scitbx::matrix::matrix_transposed_vector(
            6, 6, unit_cell.u_star_to_u_cart_linear_map().begin(),
            (u_cart_grads[grad_indices[i]]*k_ij).begin(),
            u_star_grad_ij.begin());
          for (int j=0; j < i_seqs.size(); j++) {
            scitbx::sym_mat3<double> &u_star_grad =
              (i==j ? u_star_grad_ii : u_star_grad_ij);
            for (int k=0; k < 6; k++) {
              linearised_eqns.design_matrix(first_row_i+j, ids.u_aniso+k) =
                (k > 2 ? 2*u_star_grad[k] : u_star_grad[k]);
            }
          }
        }
        else {
          CCTBX_ASSERT(ids.u_iso != -1);
          double u_iso_ii = iso_grads[grad_indices[i]]*k_ii,
            u_iso_ij = iso_grads[grad_indices[i]]*k_ij;
          for (int j=0; j < i_seqs.size(); j++) {
            linearised_eqns.design_matrix(first_row_i+j, ids.u_iso) =
              (i == j ? u_iso_ii : u_iso_ij);
          }
        }
        linearised_eqns.weights[first_row_i+i] = weight;
        linearised_eqns.deltas[first_row_i+i] = deltas_[i];
      }
    }

    double mean_u_volume;
  protected:
    af::shared<scitbx::sym_mat3<double> > u_cart_grads;
    af::shared<double> iso_grads;

    static double r3diff_to_vol(double r3diff) {
      return 4*scitbx::constants::pi*r3diff/3;
    }

    static scitbx::sym_mat3<double> calc_grad(
      adptbx::eigensystem<double> const& es)
    {
      const vec3<double> &v = es.values();
      vec3<double> vp(v[1]*v[2], v[0]*v[2], v[0]*v[1]);
      double coeff = (4*scitbx::constants::pi/3)/(2*std::sqrt(v[0]*v[1]*v[2]));
      af::tiny<double, 6> u_cart_grad;
      u_cart_grad[0] = coeff*(
        vp[0]*es.vectors(0)[0]*es.vectors(0)[0] +
        vp[1]*es.vectors(1)[0]*es.vectors(1)[0] +
        vp[2]*es.vectors(2)[0]*es.vectors(2)[0]);
      u_cart_grad[1] = coeff*(
        vp[0]*es.vectors(0)[1]*es.vectors(0)[1] +
        vp[1]*es.vectors(1)[1]*es.vectors(1)[1] +
        vp[2]*es.vectors(2)[1]*es.vectors(2)[1]);
      u_cart_grad[2] = coeff*(
        vp[0]*es.vectors(0)[2]*es.vectors(0)[2] +
        vp[1]*es.vectors(1)[2]*es.vectors(1)[2] +
        vp[2]*es.vectors(2)[2]*es.vectors(2)[2]);
      u_cart_grad[3] = 2*coeff*(
        vp[0]*es.vectors(0)[0]*es.vectors(0)[1] +
        vp[1]*es.vectors(1)[1]*es.vectors(1)[0] +
        vp[2]*es.vectors(2)[1]*es.vectors(2)[0]);
      u_cart_grad[4] = 2*coeff*(
        vp[0]*es.vectors(0)[0]*es.vectors(0)[2] +
        vp[1]*es.vectors(1)[2]*es.vectors(1)[0] +
        vp[2]*es.vectors(2)[2]*es.vectors(2)[0]);
      u_cart_grad[5] = 2*coeff*(
        vp[0]*es.vectors(0)[1]*es.vectors(0)[2] +
        vp[1]*es.vectors(1)[1]*es.vectors(1)[2] +
        vp[2]*es.vectors(2)[1]*es.vectors(2)[2]);

      return u_cart_grad;
    }
    // u_iso and u_star indices
    af::shared<std::size_t> grad_indices;
  };

}} // cctbx::adp_restraints

#endif
