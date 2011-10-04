#ifndef CCTBX_ADP_RESTRAINTS_FIXED_U_EQ_ADP_H
#define CCTBX_ADP_RESTRAINTS_FIXED_U_EQ_ADP_H

#include <cctbx/adp_restraints/adp_restraints.h>

namespace cctbx { namespace adp_restraints {

  struct fixed_u_eq_adp_proxy : public adp_restraint_proxy<1> {
    fixed_u_eq_adp_proxy() {}

    fixed_u_eq_adp_proxy(af::tiny<unsigned, 1> const &i_seq,
      double weight, double u_eq_ideal)
    : adp_restraint_proxy<1>(i_seq, weight),
      u_eq_ideal(u_eq_ideal)
    {}

    double u_eq_ideal;
  };

  class fixed_u_eq_adp : public adp_restraint_base<1> {
  public:
    fixed_u_eq_adp(scitbx::sym_mat3<double>const &u_cart,
      double weight, double u_eq_ideal_)
      : adp_restraint_base<1>(af::tiny<bool, 1>(true), weight),
      u_eq_ideal(u_eq_ideal_)
    {
      init_deltas(u_cart);
    }

    fixed_u_eq_adp(double u_iso, double weight, double u_eq_ideal_)
    : adp_restraint_base<1>(af::tiny<bool, 1>(true), weight),
      u_eq_ideal(u_eq_ideal_)
    {
      init_deltas(u_iso);
    }

    fixed_u_eq_adp(
      adp_restraint_params<double> const &params,
      fixed_u_eq_adp_proxy const &proxy)
    : adp_restraint_base<1>(params, proxy),
      u_eq_ideal(proxy.u_eq_ideal)
    {
      if (use_u_aniso[0]) {
        CCTBX_ASSERT(proxy.i_seqs[0] < params.u_cart.size());
        init_deltas(params.u_cart[proxy.i_seqs[0]]);
      }
      else {
        CCTBX_ASSERT(proxy.i_seqs[0] < params.u_iso.size());
        init_deltas(params.u_iso[proxy.i_seqs[0]]);
      }
    }

    void linearise(
      uctbx::unit_cell const &unit_cell,
      cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
      cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
      af::tiny<unsigned, 1> const &i_seqs)
    {
      linearise_1<fixed_u_eq_adp>(
        unit_cell, linearised_eqns, parameter_map, i_seqs[0], use_u_aniso[0], weight, deltas_);
    }

    void add_gradients(
      af::ref<scitbx::sym_mat3<double> > const& gradients_aniso_cart,
      af::ref<double> const& gradients_iso,
      af::tiny<unsigned, 1> const& i_seqs) const
    {
      gradients_aniso_cart[i_seqs[0]] += gradients();
    }

    static double grad_u_iso(int i) { return 1; }

    static const double* cart_grad_row(int i) {
      static const double grads_u_cart[6][6] = {
        {1./3, 1./3, 1./3, 0, 0, 0},
        {1./3, 1./3, 1./3, 0, 0, 0},
        {1./3, 1./3, 1./3, 0, 0, 0},
        {   0,    0,    0, 0, 0, 0},
        {   0,    0,    0, 0, 0, 0},
        {   0,    0,    0, 0, 0, 0}
      };
      return &grads_u_cart[i][0];
    }

    double u_eq_ideal;
  protected:

    void init_deltas(scitbx::sym_mat3<double>const &u_cart) {
      double u_eq_diff = (u_cart.trace()/3 - u_eq_ideal);
      for (int i=0; i<6; i++)
        deltas_[i] = (i < 3 ? u_eq_diff : 0);
    }

    void init_deltas(double u_iso) {
      deltas_[0] = u_iso-u_eq_ideal;
      for (int i=1; i<6; i++) deltas_[i] = 0;
    }
  };

}} // cctbx::adp_restraints

#endif
