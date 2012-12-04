#include <mmtbx/error.h>
#include <cctbx/geometry_restraints/bond.h>
#include <cctbx/geometry/geometry.h>

#include <cmath>
#include <set>
#include <iostream>

namespace mmtbx { namespace geometry_restraints {
  namespace af = scitbx::af;

  struct reference_coordinate_proxy
  {
    //! Support for shared_proxy_select.
    typedef af::tiny<unsigned, 1> i_seqs_type;

    // default initializer
    reference_coordinate_proxy () {}

    reference_coordinate_proxy(
      i_seqs_type const& i_seqs_,
      scitbx::vec3<double> ref_sites_,
      double weight_,
      double limit_=-1.0,
      bool top_out_=false)
    :
      i_seqs(i_seqs_),
      ref_sites(ref_sites_),
      weight(weight_),
      limit(limit_),
      top_out(top_out_)
    {
      if (top_out) {
        MMTBX_ASSERT(limit >= 0.0);
      }
    }

    // Support for proxy_select (and similar operations)
    reference_coordinate_proxy(
      i_seqs_type const& i_seqs_,
      reference_coordinate_proxy const& proxy)
    :
      i_seqs(i_seqs_),
      ref_sites(proxy.ref_sites),
      weight(proxy.weight),
      limit(proxy.limit),
      top_out(proxy.top_out)
    {
      if (top_out) {
        MMTBX_ASSERT(limit >= 0.0);
      }
    }

    i_seqs_type i_seqs;
    scitbx::vec3<double> ref_sites;
    double weight;
    double limit;
    bool top_out;
  };

  inline
  double
  reference_coordinate_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<reference_coordinate_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array)
  {
    double residual_sum = 0, weight, residual, limit, top;
    scitbx::vec3<double> site, ref_site, delta;
    scitbx::vec3<double> gradient;
    for (std::size_t i = 0; i < proxies.size(); i++) {
      reference_coordinate_proxy proxy = proxies[i];
      af::tiny<unsigned, 1> const& i_seqs = proxy.i_seqs;
      MMTBX_ASSERT(i_seqs[0] < sites_cart.size());
      site = sites_cart[ i_seqs[0] ];
      ref_site = proxy.ref_sites;
      weight = proxy.weight;
      delta[0] = site[0] - ref_site[0];
      delta[1] = site[1] - ref_site[1];
      delta[2] = site[2] - ref_site[2];
      if (proxy.top_out && proxy.limit >= 0.0) {
        limit = proxy.limit;
        top = weight * limit * limit;
        residual = 0.0;
        residual += top * (1.0-std::exp(-weight*delta[0]*delta[0]/top));
        residual += top * (1.0-std::exp(-weight*delta[1]*delta[1]/top));
        residual += top * (1.0-std::exp(-weight*delta[2]*delta[2]/top));
        //residual_local = ( (delta[0]*delta[0]*weight)+
        //                   (delta[1]*delta[1]*weight)+
        //                   (delta[2]*delta[2]*weight) );
        //residual = top * (1.0-std::exp(-residual_local/top));
        gradient[0] = (delta[0]*2.0*weight)
                      * std::exp(-(weight*delta[0]*delta[0])/top);
        gradient[1] = (delta[1]*2.0*weight)
                      * std::exp(-(weight*delta[1]*delta[1])/top);
        gradient[2] = (delta[2]*2.0*weight)
                      * std::exp(-(weight*delta[2]*delta[2])/top);
      }
      else {
        residual = ( (delta[0]*delta[0]*weight)+
                     (delta[1]*delta[1]*weight)+
                     (delta[2]*delta[2]*weight) );
        gradient[0] = delta[0]*2.0*weight;
        gradient[1] = delta[1]*2.0*weight;
        gradient[2] = delta[2]*2.0*weight;
      }
      residual_sum += residual;
      gradient_array[ i_seqs[0] ] += gradient;
    }
    return residual_sum;
  }
}}
