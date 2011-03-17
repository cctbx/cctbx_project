#ifndef MMTBX_GEOMETRY_RESTRAINTS_H_BOND_H
#define MMTBX_GEOMETRY_RESTRAINTS_H_BOND_H

#define IMP_SCALE 400
// sqrt(4/6)
#define IMP_SIGMA_BASE 0.81649658092772603
// pi / 180
#define PI_180 0.017453292519943295

#include <mmtbx/error.h>
#include <cctbx/geometry_restraints/bond.h>
#include <cctbx/geometry/geometry.h>

#include <cmath>
#include <set>
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
    {
      MMTBX_ASSERT((distance_cut <= 0) || (distance_cut > distance_ideal));
    }

    double distance_cut;
  };

  struct h_bond_lj_proxy
  {
    h_bond_lj_proxy() {}

    h_bond_lj_proxy(
      af::tiny<unsigned, 2> const& i_seqs_,
      double distance_ideal_,
      double distance_cut_)
    :
      i_seqs(i_seqs_),
      distance_ideal(distance_ideal_),
      distance_cut(distance_cut_)
    {
      MMTBX_ASSERT((distance_cut <= 0) || (distance_cut > distance_ideal));
    }

    af::tiny<unsigned, 2> i_seqs;
    double distance_ideal;
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
      double residual = restraint.residual();
      double grad_factor = hbond_weight;
      residual_sum += residual;
      if (gradient_array.size() != 0) {
        af::tiny<scitbx::vec3<double>, 2> gradients = restraint.gradients();
        gradient_array[ i_seqs[0] ] += gradients[0] * grad_factor;
        gradient_array[ i_seqs[1] ] += gradients[1] * grad_factor;
      }
    }
    return residual_sum;
  }


  // Switching function for Lennard-Jones-like potentials
  // I couldn't figure out how the version in the X-PLOR manual is supposed to
  // work, but NAMD has something similar:
  //   http://www.ks.uiuc.edu/Research/namd/1.5/ug/node52.html
  // which is implemented here
  inline
  double
  switch_fn (
    double R_ij,
    double R_on,
    double R_off)
  {
    if ((R_on <= 0) || (R_off <= 0)) {
      return 1.0;
    } else if (R_ij < R_on) {
      return 1.0;
    } else if (R_ij > R_off) {
      return 0.0;
    }
    double R_off_sq = R_off * R_off;
    double R_on_sq = R_on * R_on;
    double R_ij_sq = R_ij * R_ij;
    double sw = (R_off_sq - R_ij_sq) * (R_off_sq - R_ij_sq) * \
                (R_off_sq + (2 * R_ij_sq) - (3 * R_on_sq)) /
                std::pow((R_off_sq - R_on_sq), 3);
    return sw;
  }

  // from Mathematica
  //  f[r_] := ((a^2 - r^2)^2) * (a^2 + (2 * r^2) - (3 * b^2)) /((a^2 - b^2)^3)
  inline
  double
  d_switch_d_distance (
    double R_ij,
    double R_on,
    double R_off)
  {
    if ((R_on <= 0) || (R_off <= 0)) {
      return 0.0;
    } else if ((R_ij < R_on) || (R_ij > R_off)) {
      return 0.0;
    }
    double R_off_sq = R_off * R_off;
    double R_on_sq = R_on * R_on;
    double R_ij_sq = R_ij * R_ij;
    double R_off_sq_minus_R_on_sq_cub = std::pow((R_off_sq - R_on_sq), 3);
    double d_sw_d_R = ((4*R_ij*(R_off_sq-R_ij_sq)*(R_off_sq-R_ij_sq)) /
                        R_off_sq_minus_R_on_sq_cub) -
         ((4*R_ij*(R_off_sq-R_ij_sq)*(R_off_sq-(3*R_on_sq)+(2*R_ij_sq))) /
                        R_off_sq_minus_R_on_sq_cub);
    return d_sw_d_R;
  }

  inline
  double
  eval_lennard_jones_energy(
    double r_ij,
    double sigma,
    int a,
    int b)
  {
    MMTBX_ASSERT(r_ij > 0);
    return (std::pow(sigma / r_ij, a) - std::pow(sigma / r_ij, b));
  }

  inline
  double
  eval_lennard_jones_deriv(
    double r_ij,
    double sigma,
    int a,
    int b)
  {
    MMTBX_ASSERT(r_ij > 0);
    return (-a * std::pow(sigma, a) / std::pow(r_ij, a+1)) +
            (b * std::pow(sigma, b) / std::pow(r_ij, b+1));
  }

  inline
  double
  h_bond_lennard_jones_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<h_bond_lj_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    double falloff_distance=0.05,
    double sigma_base=0.81649658092772603,
    int a=6,
    int b=4,
    double scale=1.0,
    bool use_finite_differences=false)
  {
    using namespace cctbx::geometry;
    double sum = 0.0;
    for (std::size_t i = 0; i < proxies.size(); i++) {
      h_bond_lj_proxy proxy = proxies[i];
      af::tiny<unsigned, 2> const& i_seqs = proxy.i_seqs;
      af::tiny<scitbx::vec3<double>, 2> bond_sites;
      bond_sites[0] = sites_cart[i_seqs[0]];
      bond_sites[1] = sites_cart[i_seqs[1]];
      distance<double> r_da(bond_sites);
      double r_ij = r_da.distance_model;
      double sigma = sigma_base * proxy.distance_ideal;
      double weight = scale; // * hbond_weight;
      double residual = weight * eval_lennard_jones_energy(r_ij, sigma, a, b);
      double sw = 1;
      double R_off = proxy.distance_cut + falloff_distance;
      if (proxy.distance_cut > 0) {
        sw = switch_fn(r_ij, proxy.distance_cut, R_off);
      }
      sum += residual * sw;
      if (gradient_array.size() != 0) {
        if (use_finite_differences) {
          for (unsigned j = 0; j < 2; j++) {
            for (unsigned k = 0; k < 3; k++) {
              bond_sites[j][k] -= 0.0001;
              distance<double> r_da_1(bond_sites);
              double r_ij_1 = r_da_1.distance_model;
              double e1 = eval_lennard_jones_energy(r_ij_1, sigma, a, b) *
                          switch_fn(r_ij_1, proxy.distance_cut, R_off);
              bond_sites[j][k] += 0.0002;
              distance<double> r_da_2(bond_sites);
              double r_ij_2 = r_da_2.distance_model;
              double e2 = eval_lennard_jones_energy(r_ij_2, sigma, a, b)*
                          switch_fn(r_ij_2, proxy.distance_cut, R_off);
              bond_sites[j][k] -= 0.0001;
              double grad_k = weight * (e2 - e1) / 0.0002;
              gradient_array[ i_seqs[j] ][k] += grad_k;
            }
          }
        } else {
          af::tiny<scitbx::vec3<double>, 2> d_R_d_xyz = \
            r_da.d_distance_d_sites();
          double d_lj_d_R = (-a * std::pow(sigma, a) / std::pow(r_ij, a+1)) +
                            (b * std::pow(sigma, b) / std::pow(r_ij, b+1));
          double d_sw_d_R = d_switch_d_distance(r_ij,proxy.distance_cut,R_off);
          for (unsigned j = 0; j < 2; j++) {
            scitbx::vec3<double> d_lj_d_xyz = weight * d_lj_d_R * d_R_d_xyz[j];
            scitbx::vec3<double> d_sw_d_xyz = d_sw_d_R * d_R_d_xyz[j];
            gradient_array[ i_seqs[j] ] -= (sw * d_lj_d_xyz) +
                                           (residual * d_sw_d_xyz);
          }
        }
      }
    }
    return sum;
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
    double R_off = distance_cut + falloff_distance;
    if (bond_dist > R_off) {
      return 0.0;
    }
    double sigma_over_dist = distance_ideal * IMP_SIGMA_BASE / bond_dist;
    double s_over_d_sq = sigma_over_dist * sigma_over_dist;
    double residual = weight * ((s_over_d_sq * s_over_d_sq * s_over_d_sq) - \
                                (s_over_d_sq * s_over_d_sq)) * \
                      std::pow(std::cos(delta_theta * PI_180), 4);
    residual *= switch_fn(bond_dist, distance_cut, R_off);
    return residual;
  }

  inline
  double
  h_bond_implicit_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<h_bond_implicit_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    double falloff_distance=0.05)
  {
    using namespace cctbx::geometry;
    double residual_sum = 0.0;
    for (std::size_t i = 0; i < proxies.size(); i++) {
      h_bond_implicit_proxy proxy = proxies[i];
      af::tiny<scitbx::vec3<double>, 3> sites;
      af::tiny<unsigned, 3> const& i_seqs = proxy.i_seqs;
      for (unsigned j = 0; j < 3; j++) {
        sites[j] = sites_cart[i_seqs[j]];
      }
      double weight = proxy.weight * IMP_SCALE;
      angle<double> theta(sites);
      if (!theta.have_angle_model) continue;
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
        af::tiny<scitbx::vec3<double>, 2> bond_sites;
        bond_sites[0] = sites[0];
        bond_sites[1] = sites[1];
        distance<double> r_da(bond_sites);
        af::tiny<scitbx::vec3<double>, 2> d_R_d_xyz = \
          r_da.d_distance_d_sites();
        af::tiny<scitbx::vec3<double>, 3> d_theta_d_xyz = \
          theta.d_angle_d_sites();
        // f(R) = L-J term
        // g(theta) = angle term
        // h(R) = switching function
        double R = r_da.distance_model;
        double R_on = proxy.distance_cut;
        double R_off = R_on + falloff_distance;
        double h_R = switch_fn(R, R_on, R_off);
        double d_h_d_R = d_switch_d_distance(R, R_on, R_off) / weight;
        double sigma = proxy.distance_ideal * IMP_SIGMA_BASE;
        double f_R = weight * (std::pow(sigma/R, 6) - std::pow(sigma/R, 4));
        double d_f_d_R = (-6 * std::pow(sigma, 6) / std::pow(R, 7)) +
                         (4 * std::pow(sigma, 4) / std::pow(R, 5));
        double cos_delta_theta = std::cos(delta_theta * PI_180);
        double g_theta = weight * std::pow(cos_delta_theta, 4);
        double d_g_d_theta = 4 * std::pow(cos_delta_theta, 3) *
                             std::sin(delta_theta * PI_180);
        for (unsigned j = 0; j < 3; j++) {
          scitbx::vec3<double> d_g_d_xyz = d_g_d_theta * d_theta_d_xyz[j];
          scitbx::vec3<double> grads = f_R * h_R * d_g_d_xyz * PI_180;
          if (j != 2) {
            scitbx::vec3<double> d_h_d_xyz = d_h_d_R * d_R_d_xyz[j];
            grads -= f_R * g_theta * d_h_d_xyz;
            scitbx::vec3<double> d_f_d_xyz = d_f_d_R * d_R_d_xyz[j];
            grads -= h_R * g_theta * d_f_d_xyz; // XXX why minus?
          }
          gradient_array[i_seqs[j]] += grads;
        }
      }
    }
    return residual_sum;
  }

  af::shared<std::set<unsigned> > simple_hbonds_as_simple_bonds (
    af::const_ref<h_bond_simple_proxy> const& proxies)
  {
    af::shared<std::set<unsigned> > bonded_pairs;
    for (unsigned i = 0; i < proxies.size(); i++) {
      std::set<unsigned> pair;
      pair.insert(proxies[i].i_seqs[0]);
      pair.insert(proxies[i].i_seqs[1]);
      bonded_pairs.push_back(pair);
    }
    return bonded_pairs;
  }

  af::shared<std::set<unsigned> > lj_hbonds_as_simple_bonds (
    af::const_ref<h_bond_lj_proxy> const& proxies)
  {
    af::shared<std::set<unsigned> > bonded_pairs;
    for (unsigned i = 0; i < proxies.size(); i++) {
      std::set<unsigned> pair;
      pair.insert(proxies[i].i_seqs[0]);
      pair.insert(proxies[i].i_seqs[1]);
      bonded_pairs.push_back(pair);
    }
    return bonded_pairs;
  }

  af::shared<std::set<unsigned> > implicit_hbonds_as_simple_bonds (
    af::const_ref<h_bond_implicit_proxy> const& proxies)
  {
    af::shared<std::set<unsigned> > bonded_pairs;
    for (unsigned i = 0; i < proxies.size(); i++) {
      std::set<unsigned> pair;
      pair.insert(proxies[i].i_seqs[0]);
      pair.insert(proxies[i].i_seqs[1]);
      bonded_pairs.push_back(pair);
    }
    return bonded_pairs;
  }

}} // namespace mmtbx::geometry_restraints

#endif
