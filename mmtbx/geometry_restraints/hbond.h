#ifndef MMTBX_GEOMETRY_RESTRAINTS_H_BOND_H
#define MMTBX_GEOMETRY_RESTRAINTS_H_BOND_H

// sqrt(2/3)
#define SIGMA_BASE 0.81649658092772603
// 2 * pi / 180
#define TWOPI_180 0.034906585039886591

#include <mmtbx/error.h>
#include <cctbx/geometry_restraints/bond.h>
#include <cctbx/geometry/geometry.h>

#include <cmath>
#include <iostream>

namespace mmtbx { namespace geometry_restraints {

  namespace af = scitbx::af;
  using cctbx::geometry_restraints::bond;
  using cctbx::geometry_restraints::bond_simple_proxy;

  // H:O or N:O pseudo-bond (uses cctbx bond energy internally)
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

  // N:O=C with both distance and angle terms (see below)
  struct h_bond_implicit_proxy
  {
    h_bond_implicit_proxy() {}

    h_bond_implicit_proxy(
      af::tiny<unsigned, 3> const& i_seqs_,
      double distance_ideal_,
      double distance_cut_,
      double theta_low_,
      double theta_high_,
      double weight_)
    :
      i_seqs(i_seqs_),
      distance_ideal(distance_ideal_),
      distance_cut(distance_cut_),
      theta_low(theta_low_),
      theta_high(theta_high_),
      weight(weight_)
    {
      MMTBX_ASSERT((weight>=0) && (distance_ideal>=0) && (distance_cut>=0));
    }

    af::tiny<unsigned, 3> i_seqs;
    double distance_ideal;
    double distance_cut;
    double theta_low;
    double theta_high;
    double weight;
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

  // Fabiola et al. (2002) Protein Sci. 11:1415-23
  // http://www.ncbi.nlm.nih.gov/pubmed/12021440
  inline
  double
  residual_implicit(
    af::tiny<scitbx::vec3<double>, 3> sites,
    double distance_ideal,
    double distance_cut,
    double weight,
    double delta_theta,
    double falloff_distance)
  {
    double bond_dist = (sites[1] - sites[0]).length();
    //std::cout << bond_dist << "\n";
    double delta_cut = bond_dist - distance_cut;
    if (delta_cut > falloff_distance) {
      return 0.0;
    }
    double sigma_over_dist = distance_ideal * SIGMA_BASE / bond_dist;
    double s_over_d_sq = sigma_over_dist * sigma_over_dist;
    double residual = weight * ((s_over_d_sq * s_over_d_sq * s_over_d_sq) - \
                                (s_over_d_sq * s_over_d_sq)) * \
                      std::pow(std::cos(delta_theta * TWOPI_180), 4);
    if (delta_cut > 0) {
      residual *= (falloff_distance - delta_cut) / falloff_distance;
    }
    return residual;
  }

  // gradients by finite differences
  inline
  double
  h_bond_implicit_residual_sum_fd(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<h_bond_implicit_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    double falloff_distance=0.05,
    double epsilon=0.0001)
  {
    using cctbx::geometry::angle;
    double residual_sum = 0.0;
    double two_epsilon = epsilon * 2;
    for (std::size_t i = 0; i < proxies.size(); i++) {
      h_bond_implicit_proxy proxy = proxies[i];
      af::tiny<scitbx::vec3<double>, 3> sites;
      af::tiny<unsigned, 3> const& i_seqs = proxy.i_seqs;
      for (unsigned j = 0; j < 3; j++) {
        sites[j] = sites_cart[i_seqs[j]];
      }
      double weight = proxy.weight;
      angle<double> theta(sites);
      double angle_ideal, delta_theta;
      double delta_high = theta.angle_model - proxy.theta_high;
      double delta_low = theta.angle_model - proxy.theta_low;
      if ((proxy.theta_low < 0) || (fabs(delta_high) < fabs(delta_low))) {
        delta_theta = delta_high;
        angle_ideal = proxy.theta_high;
      } else {
        delta_theta = delta_low;
        angle_ideal = proxy.theta_low;
      }
      residual_sum += residual_implicit(
        sites,
        proxy.distance_ideal,
        proxy.distance_cut,
        weight,
        delta_theta,
        falloff_distance);
      if (gradient_array.size() != 0) {
        for (unsigned j = 0; j < 3; j++) {
          for (unsigned k = 0; k < 3; k++) {
            sites[j][k] -= epsilon;
            angle<double> theta_1(sites);
            double delta_theta_1 = theta_1.angle_model - angle_ideal;
            double e1 = residual_implicit(
              sites,
              proxy.distance_ideal,
              proxy.distance_cut,
              weight,
              delta_theta_1,
              falloff_distance);
            sites[j][k] += two_epsilon;
            angle<double> theta_2(sites);
            double delta_theta_2 = theta_2.angle_model - angle_ideal;
            double e2 = residual_implicit(
              sites,
              proxy.distance_ideal,
              proxy.distance_cut,
              weight,
              delta_theta_2,
              falloff_distance);
            gradient_array[i_seqs[j]][k] += (e2 - e1) / two_epsilon;
            sites[j][k] -= epsilon;
          }
        }
      }
    }
    return residual_sum;
  }

}} // namespace mmtbx::geometry_restraints

#endif
