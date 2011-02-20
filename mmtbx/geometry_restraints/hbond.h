#ifndef MMTBX_GEOMETRY_RESTRAINTS_H_BOND_H
#define MMTBX_GEOMETRY_RESTRAINTS_H_BOND_H

#include <cctbx/geometry_restraints/bond.h>

namespace mmtbx { namespace geometry_restraints {

  namespace af = scitbx::af;
  using cctbx::geometry_restraints::bond;
  using cctbx::geometry_restraints::bond_simple_proxy;

  struct h_bond_simple_proxy : bond_simple_proxy
  {
    h_bond_simple_proxy() {}

    h_bond_simple_proxy(
      af::tiny<unsigned, 2> const& i_seqs_,
      double distance_ideal_,
      double distance_cut_,
      double weight_,
      double slack_=0)
    :
      bond_simple_proxy(i_seqs_, distance_ideal_, weight_, slack_),
      distance_cut(distance_cut_)
    {}

    double distance_cut;
  };

  inline
  double
  h_bond_simple_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<h_bond_simple_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    double hbond_weight=1.0,
    double falloff_distance=0.05)
  {
    double residual_sum = 0;
    for (std::size_t i = 0; i < proxies.size(); i++) {
      h_bond_simple_proxy proxy = proxies[i];
      af::tiny<scitbx::vec3<double>, 2> sites;
      af::tiny<unsigned, 2> const& i_seqs = proxy.i_seqs;
      sites[0] = sites_cart[ i_seqs[0] ];
      sites[1] = sites_cart[ i_seqs[1] ];
      bond restraint(sites, proxy.distance_ideal, proxy.weight, proxy.slack);
      double cutoff_weight = hbond_weight;
      double cutoff_delta = restraint.distance_model - proxy.distance_cut;
      if (proxy.distance_cut <= 0) {
        cutoff_delta = 1.0;
      } else if (cutoff_delta > falloff_distance) {
        continue;
      } else if (cutoff_delta > 0) {
        cutoff_weight *= (falloff_distance - cutoff_delta) / falloff_distance;
      }
      double residual = restraint.residual() * cutoff_weight;
      residual_sum += residual;
      if (gradient_array.size() != 0) {
        af::tiny<scitbx::vec3<double>, 2> gradients = restraint.gradients();
        gradient_array[ i_seqs[0] ] += gradients[0] * cutoff_weight;
        gradient_array[ i_seqs[1] ] += gradients[1] * cutoff_weight;
      }
    }
    return residual_sum;
  }

}} // namespace mmtbx::geometry_restraints

#endif
