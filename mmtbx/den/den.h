#include <mmtbx/error.h>
#include <cctbx/geometry_restraints/bond.h>
#include <cctbx/geometry/geometry.h>

#include <cmath>
#include <set>
#include <iostream>

namespace mmtbx { namespace den {

  namespace af = scitbx::af;
  using cctbx::geometry_restraints::bond;
  using cctbx::geometry_restraints::bond_simple_proxy;

  struct den_simple_proxy : bond_simple_proxy
  {
    typedef af::tiny<unsigned, 2> i_seqs_type;

    den_simple_proxy() {}

    den_simple_proxy(
      i_seqs_type const& i_seqs_,
      double eq_distance_,
      double eq_distance_start_,
      double weight_)
    :
      i_seqs(i_seqs_),
      eq_distance(eq_distance_),
      eq_distance_start(eq_distance_start_),
      weight(weight_)
    {
      MMTBX_ASSERT((eq_distance > 0) && (eq_distance_start > 0));
    }

    // Support for proxy_select (and similar operations)
    den_simple_proxy(
      i_seqs_type const& i_seqs_,
      den_simple_proxy const& proxy)
    :
      i_seqs(i_seqs_),
      eq_distance(proxy.eq_distance),
      eq_distance_start(proxy.eq_distance_start),
      weight(proxy.weight)
    {
      MMTBX_ASSERT((eq_distance > 0) && (eq_distance_start > 0));
    }

    i_seqs_type i_seqs;
    double eq_distance;
    double eq_distance_start;
    double weight;
  };

  inline
  double
  den_simple_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<den_simple_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    double den_weight=1.0)
  {
    double residual_sum = 0;
    double slack = 0.0;
    unsigned n_sites = sites_cart.size();
    for (std::size_t i = 0; i < proxies.size(); i++) {
      den_simple_proxy proxy = proxies[i];
      af::tiny<scitbx::vec3<double>, 2> sites;
      af::tiny<unsigned, 2> const& i_seqs = proxy.i_seqs;
      sites[0] = sites_cart[ i_seqs[0] ];
      sites[1] = sites_cart[ i_seqs[1] ];
      MMTBX_ASSERT((i_seqs[0] < n_sites) && (i_seqs[1] < n_sites));
      //bond restraint(sites, proxy.eq_distance, proxy.weight, slack);
      bond restraint(sites, proxy.eq_distance, den_weight, slack);
      double residual = restraint.residual();
      //double grad_factor = den_weight;
      residual_sum += residual;
      if (gradient_array.size() != 0) {
        af::tiny<scitbx::vec3<double>, 2> gradients = restraint.gradients();
        //weight is now handled at the residual level
        gradient_array[ i_seqs[0] ] += gradients[0];// * grad_factor;
        gradient_array[ i_seqs[1] ] += gradients[1];// * grad_factor;
      }
    }
    return residual_sum;
  }

  inline
  void
  den_update_eq_distances(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::ref<den_simple_proxy> const& proxies,
    double gamma,
    double kappa)
  {
    double slack = 0.0;
    for (std::size_t i = 0; i < proxies.size(); i++) {
      den_simple_proxy proxy = proxies[i];
      af::tiny<scitbx::vec3<double>, 2> sites;
      af::tiny<unsigned, 2> const& i_seqs = proxy.i_seqs;
      sites[0] = sites_cart[ i_seqs[0] ];
      sites[1] = sites_cart[ i_seqs[1] ];
      bond restraint(sites, proxy.eq_distance, proxy.weight, slack);
      double distance_model = restraint.distance_model;
      double new_eq_dist = ((1.0-kappa)*proxy.eq_distance) +
                           ( kappa *
                             ( (gamma*distance_model) +
                               (1.0-gamma)*proxy.eq_distance_start ));
      proxies[i].eq_distance = new_eq_dist;
    }
  }
}}
