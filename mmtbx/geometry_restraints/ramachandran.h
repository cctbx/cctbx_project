#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
//#include <boost/python/return_value_policy.hpp>
//#include <boost/python/return_by_value.hpp>
#include <boost/optional.hpp>

#include <mmtbx/error.h>
#include <cctbx/geometry_restraints/dihedral.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/math/interpolation.h>

#include <cmath>
#include <iostream>

namespace mmtbx { namespace geometry_restraints {
  using cctbx::geometry_restraints::dihedral;
  using cctbx::geometry_restraints::dihedral_proxy;
  namespace af = scitbx::af;

  struct phi_psi_proxy
  {
    typedef af::tiny<unsigned, 5> i_seqs_type;

    i_seqs_type i_seqs;
    std::string residue_name;
    std::string residue_type;
    size_t residue_index;

    // default initializer
    phi_psi_proxy () {}

    phi_psi_proxy(
      i_seqs_type i_seqs_,
      std::string const& residue_name_,
      std::string const& residue_type_,
      size_t residue_index_=1)
    :
      i_seqs(i_seqs_),
      residue_name(residue_name_),
      residue_type(residue_type_),
      residue_index(residue_index_)
    {
      MMTBX_ASSERT(residue_index > 0);
    }

    phi_psi_proxy(
      i_seqs_type const& i_seqs_,
      phi_psi_proxy const& proxy)
    :
      i_seqs(i_seqs_),
      residue_name(proxy.residue_name),
      residue_type(proxy.residue_type),
      residue_index(proxy.residue_index)
    {
      MMTBX_ASSERT(residue_index > 0);
    }

    af::shared<unsigned> get_i_seqs() {
      af::shared<unsigned> result;
      for (int i=0; i<5; i++) result.push_back(i_seqs[i]);
      return result;
    }

  };

  namespace gr = cctbx::geometry_restraints;

  double convert_angle (double theta) {
    if (theta > 180) {
      theta = -360 + theta;
    } else if (theta < -180) {
      theta = 360 + theta;
    }
    return theta;
  }

  // COOT-like restraints (rama_potential=emsley)
  class lookup_table
  {
    public:
      af::versa<double, af::flex_grid<> > plot;
      double values_max;

      lookup_table (
        af::const_ref< double > values,
        int n_angles,
        double scale_allowed=1.0)
      {
        MMTBX_ASSERT(values.size() == (n_angles * n_angles));
        MMTBX_ASSERT(scale_allowed > 0.0);
        af::flex_grid<>::index_type fg_origin;
        af::flex_grid<>::index_type fg_last;
        for (unsigned i = 0; i < 2; i++) {
          fg_origin.push_back(0);
          fg_last.push_back(n_angles);
        }
        plot = af::versa<double, af::flex_grid<> >(
          af::flex_grid<>(fg_origin, fg_last, true));
        values_max = 0.0;
        for (unsigned i = 0; i < values.size(); i++) {
          if (values[i] > 0.0) {
            plot[i] = values[i] * scale_allowed;
          } else {
            plot[i] = values[i];
          }
          if (plot[i] > values_max) {
            values_max = plot[i];
          }
        }
      }

      double get_score (
        double phi,
        double psi,
        bool use_splines=false)
      {
        using scitbx::math::interpolate_at_point;
        phi = convert_angle(phi);
        psi = convert_angle(psi);
        //std::cout << phi << " " << psi << "\n";
        MMTBX_ASSERT((phi <= 180.0) && (phi >= -180.0));
        MMTBX_ASSERT((psi <= 180.0) && (psi >= -180.0));
        int phi_1 = (int) floor(phi);
        int phi_2 = (int) ceil(phi);
        int psi_1 = (int) floor(psi);
        int psi_2 = (int) ceil(psi);
        if ((phi_1 % 2) == 0) {
          if (phi_2 == phi_1)
            phi_2 += 1;
          phi_1 -= 1;
        } else if ((phi_2 % 2) == 0) {
          phi_2 += 1;
        }
        if ((psi_1 % 2) == 0) {
          if (psi_2 == psi_1)
            psi_2 += 1;
          psi_1 -= 1;
        } else if ((psi_2 % 2) == 0) {
          psi_2 += 1;
        }
        double r_phi_psi = 0;
        if (use_splines) {
          double x_phi = (phi - (double) phi_1) / 2.0;
          double x_psi = (psi - (double) psi_1) / 2.0;
          af::tiny<double, 4> phi_n(0.0);
          for (int i = -1; i < 3; i++) {
            int phi_i = phi_1 + (i * 2);
            af::tiny<double, 4> psi_n(0.0);
            for (int j = -1; j < 3; j++) {
              int psi_j = psi_1 + (j * 2);
              psi_n[j+1] = get_point(phi_i, psi_j);
            }
            phi_n[i+1] = interpolate_at_point(psi_n, x_psi);
          }
          r_phi_psi = interpolate_at_point(phi_n, x_phi);
        } else if (phi_2 == phi_1) {
          if (psi_2 == psi_1) {
            r_phi_psi = get_point(phi_1, psi_1);
          } else {
            double r11 = get_point(phi_1, psi_1);
            double r12 = get_point(phi_1, psi_2);
            r_phi_psi = (((psi-psi_1)*r12)+((psi_2-psi)*r11)) / (psi_2-psi_1);
          }
        } else if (psi_2 == psi_1) {
          double r11 = get_point(phi_1, psi_1);
          double r21 = get_point(phi_2, psi_1);
          r_phi_psi = (((phi-phi_1)*r21)+((phi_2-phi)*r11)) / (phi_2-phi_1);
        } else {
          double r11 = get_point(phi_1, psi_1);
          double r12 = get_point(phi_1, psi_2);
          double r21 = get_point(phi_2, psi_1);
          double r22 = get_point(phi_2, psi_2);
          double d_phi_d_psi = (double) (phi_2 - phi_1) * (psi_2 - psi_1);
          MMTBX_ASSERT(d_phi_d_psi != 0);
          r_phi_psi = ((r11/d_phi_d_psi) * (phi_2-phi) * (psi_2-psi)) +\
                      ((r21/d_phi_d_psi) * (phi-phi_1) * (psi_2-psi)) +\
                      ((r12/d_phi_d_psi) * (phi_2-phi) * (psi-psi_1)) +\
                      ((r22/d_phi_d_psi) * (phi-phi_1) * (psi-psi_1));
        }
        return r_phi_psi;
      }

      // score inverted and adjusted to have a minimum value of 0
      double get_energy (
        double phi,
        double psi)
      {
        double score = get_score(phi, psi);
        return - (score - values_max);
      }

      double
      compute_gradients (
        af::ref<scitbx::vec3<double> > const& gradient_array,
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        phi_psi_proxy const& proxy,
        double weight=1.0,
        double epsilon=0.1)
      {
        MMTBX_ASSERT(gradient_array.size() == sites_cart.size());
        MMTBX_ASSERT(epsilon > 0.0);
        af::tiny<scitbx::vec3<double>, 4> phi_sites;
        af::tiny<scitbx::vec3<double>, 4> psi_sites;
        af::tiny<unsigned, 5> const i_seqs = proxy.i_seqs;
        for (unsigned i = 0; i < 4; i++) {
          phi_sites[i] = sites_cart[i_seqs[i]];
          psi_sites[i] = sites_cart[i_seqs[i+1]];
        }
        dihedral phi(phi_sites, 0, 1.0);
        dihedral psi(psi_sites, 0, 1.0);
        double phi_deg = phi.angle_model;
        double psi_deg = psi.angle_model;
        double residual = get_energy(phi_deg, psi_deg);
        double r_phi_1 = get_energy(phi_deg - epsilon, psi_deg);
        double r_phi_2 = get_energy(phi_deg + epsilon, psi_deg);
        double d_r_d_phi = (r_phi_2 - r_phi_1) / (epsilon * 2);
        double r_psi_1 = get_energy(phi_deg, psi_deg - epsilon);
        double r_psi_2 = get_energy(phi_deg, psi_deg + epsilon);
        double d_r_d_psi = (r_psi_2 - r_psi_1) / (epsilon * 2);
        af::tiny<scitbx::vec3<double>, 4> d_phi_d_xyz = - phi.grad_delta();
        af::tiny<scitbx::vec3<double>, 4> d_psi_d_xyz = - psi.grad_delta();
        for (unsigned k = 0; k < 5; k++) {
          std::size_t i_seq = i_seqs[k];
          if (k < 4) {
            gradient_array[i_seq] += d_r_d_phi * d_phi_d_xyz[k] * weight;
          }
          if (k > 0) {
            gradient_array[i_seq] += d_r_d_psi * d_psi_d_xyz[k-1] * weight;
          }
        }
        return residual * weight;
        // this is to get magnitudes of gradients instead of function value
        // return (d_r_d_phi*d_r_d_phi+d_r_d_psi*d_r_d_psi) * weight;
      }

    private :
      double get_point (
        int phi,
        int psi)
      {
        phi = (int) convert_angle(phi);
        psi = (int) convert_angle(psi);
        MMTBX_ASSERT((phi < 180) && (phi > -180));
        MMTBX_ASSERT((psi < 180) && (psi > -180));
        MMTBX_ASSERT((abs(phi % 2) == 1) && (abs(psi % 2) == 1));
        int i = (int) (phi + 179) / 2;
        int j = (int) (psi + 179) / 2;
        return plot(i,j);
      }
  };

  template <typename FloatType>
  af::tiny<FloatType, 3>
    target_phi_psi(af::const_ref<scitbx::vec3<double> > const& rama_table,
                   af::const_ref<scitbx::vec3<double> > const& sites_cart,
                   phi_psi_proxy const& proxy)
  {
    af::tiny<unsigned, 5> const i_seqs = proxy.i_seqs;
    af::tiny<scitbx::vec3<double>, 4> phi_sites;
    af::tiny<scitbx::vec3<double>, 4> psi_sites;
    for (unsigned i = 0; i < 4; i++) {
      phi_sites[i] = sites_cart[i_seqs[i]];
      psi_sites[i] = sites_cart[i_seqs[i+1]];
    }
    dihedral phi1(phi_sites, 0, 1.0);
    dihedral psi1(psi_sites, 0, 1.0);
    double phi_deg = phi1.angle_model;
    double psi_deg = psi1.angle_model;
    double phi_t=phi_deg;
    double psi_t=psi_deg;
    double dist_to_current = 1.e+9;
    double score_current = 0;
    af::shared<double> distances;
    distances.resize(rama_table.size(), 0);
    // This is just to find out score of the current point
    for(int i=0; i<rama_table.size(); i++) {
      scitbx::vec3<double> point = rama_table[i];
      double d1 = gr::angle_delta_deg(point[0],phi_deg);
      double d2 = gr::angle_delta_deg(point[1],psi_deg);
      double d = std::sqrt(d1 * d1 + d2 * d2);
      distances[i] = d;
      if(d<dist_to_current) {
        dist_to_current = d;
        score_current = point[2];
      }
    }
    // Now we just looking for the nearest point with better score than current
    // point. then return phi-psi and distance to that point.
    // distance is used to determine weight for gradient calculation.
    double dist_to_allowed = 1.e+9;
    for(int i=0; i<rama_table.size(); i++) {
      scitbx::vec3<double> point = rama_table[i];
      double d = distances[i];
      if(point[2] >= score_current && d<dist_to_allowed) {
        dist_to_allowed = d;
        phi_t = point[0];
        psi_t = point[1];
      }
    }
    return af::tiny<double, 3> (phi_t,psi_t,dist_to_allowed) ;
  };

  af::shared<scitbx::vec3<double> >
  phi_psi_targets(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<phi_psi_proxy> const& proxies,
    af::const_ref<scitbx::vec3<double> > const& general_table,
    af::const_ref<scitbx::vec3<double> > const& gly_table,
    af::const_ref<scitbx::vec3<double> > const& cispro_table,
    af::const_ref<scitbx::vec3<double> > const& transpro_table,
    af::const_ref<scitbx::vec3<double> > const& prepro_table,
    af::const_ref<scitbx::vec3<double> > const& ileval_table)
  {
    af::shared<scitbx::vec3<double> > result;
    result.resize(proxies.size(), scitbx::vec3<double>(0,0,0));
    for (std::size_t i=0; i<proxies.size(); i++) {
      phi_psi_proxy const& proxy = proxies[i];
      af::tiny<double, 3> r;
      if (proxy.residue_type == "general") {
        r = target_phi_psi<double>(general_table, sites_cart, proxy);}
      else if (proxy.residue_type == "glycine") {
        r = target_phi_psi<double>(gly_table, sites_cart, proxy);}
      else if (proxy.residue_type == "cis-proline") {
        r = target_phi_psi<double>(cispro_table, sites_cart, proxy);}
      else if (proxy.residue_type == "trans-proline") {
        r = target_phi_psi<double>(transpro_table, sites_cart, proxy);}
      else if (proxy.residue_type == "pre-proline") {
        r = target_phi_psi<double>(prepro_table, sites_cart, proxy);}
      else if (proxy.residue_type == "isoleucine or valine") {
        r = target_phi_psi<double>(ileval_table, sites_cart, proxy);}
      else {
        std::string error_msg;
        error_msg = "Wrong proxy_type in Ramachandran proxy: '";
        error_msg += proxy.residue_type + "'";
        throw error(error_msg);
      }
      result[i]=r;
    }
      return result;
  };

  // QUANTA-like harmonic restraints (rama_potential=oldfield)
  double
  ramachandran_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<phi_psi_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    af::ref<scitbx::vec3<double> > const& phi_psi_targets,
    af::tiny<double, 4> const& weights,
    af::ref<double > const& residuals_array)
  {
    MMTBX_ASSERT(gradient_array.size() == sites_cart.size());
    MMTBX_ASSERT(residuals_array.size() == proxies.size());
    MMTBX_ASSERT(phi_psi_targets.size() == proxies.size());
    double res_sum = 0;
    for (std::size_t i=0; i<proxies.size(); i++) {
      phi_psi_proxy const& proxy = proxies[i];
      // r: phi target, psi target, distance to allowed.
      scitbx::vec3<double> r = phi_psi_targets[i];
      double weight = weights[0];
      if(std::abs(weight) < 1e-7) {
        weight = weights[1] * std::max(weights[2], std::min(r[2], weights[3]));
      }
      af::tiny<scitbx::vec3<double>, 4> phi_sites;
      af::tiny<scitbx::vec3<double>, 4> psi_sites;
      af::tiny<unsigned, 5> const i_seqs = proxy.i_seqs;
      for (unsigned j = 0; j < 4; j++) {
        phi_sites[j] = sites_cart[i_seqs[j]];
        psi_sites[j] = sites_cart[i_seqs[j+1]];
      }
      dihedral phi(phi_sites, r[0], weight);
      dihedral psi(psi_sites, r[1], weight);
      double res = phi.residual()+psi.residual();
      res_sum += res;
      residuals_array[i] = res;
      af::tiny<scitbx::vec3<double>, 4> d_phi_d_xyz = phi.gradients();
      af::tiny<scitbx::vec3<double>, 4> d_psi_d_xyz = psi.gradients();
      gradient_array[i_seqs[0]] += d_phi_d_xyz[0];
      gradient_array[i_seqs[1]] += d_phi_d_xyz[1]+d_psi_d_xyz[0];
      gradient_array[i_seqs[2]] += d_phi_d_xyz[2]+d_psi_d_xyz[1];
      gradient_array[i_seqs[3]] += d_phi_d_xyz[3]+d_psi_d_xyz[2];
      gradient_array[i_seqs[4]] += d_psi_d_xyz[3];
    }
    return res_sum;
  }


}} // namespace mmtbx::geometry_restraints
