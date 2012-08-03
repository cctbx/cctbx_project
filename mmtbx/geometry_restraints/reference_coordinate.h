#include <mmtbx/error.h>
#include <cctbx/geometry_restraints/bond.h>
#include <cctbx/geometry/geometry.h>
//#include <scitbx/math/distance_difference.h>

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
      double weight_)
    :
      i_seqs(i_seqs_),
      ref_sites(ref_sites_),
      weight(weight_)
    {}

    // Support for proxy_select (and similar operations)
    reference_coordinate_proxy(
      i_seqs_type const& i_seqs_,
      reference_coordinate_proxy const& proxy)
    :
      i_seqs(i_seqs_),
      ref_sites(proxy.ref_sites),
      weight(proxy.weight)
    {}

    i_seqs_type i_seqs;
    scitbx::vec3<double> ref_sites;
    double weight;
  };

  inline
  double
  reference_coordinate_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<reference_coordinate_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array)
  {
    double residual_sum = 0, weight;
    //unsigned n_sites = sites_cart.size();
    //af::tiny<scitbx::vec3<double>, 2> sites;
    scitbx::vec3<double> site, ref_site, delta;
    //af::versa< FloatType, af::c_grid<2> > ddm;
    //unsigned i_seq, ref_i_seq;
    scitbx::vec3<double> gradient;
    for (std::size_t i = 0; i < proxies.size(); i++) {
      reference_coordinate_proxy proxy = proxies[i];
      af::tiny<unsigned, 1> const& i_seqs = proxy.i_seqs;
      //af::tiny<unsigned, 1> const& ref_i_seqs = proxy.ref_i_seqs;
      site = sites_cart[ i_seqs[0] ];
      ref_site = proxy.ref_sites;
      weight = proxy.weight;
      //ddm = distance_difference_matrix(site, ref_site);
      delta[0] = site[0] - ref_site[0];
      delta[1] = site[1] - ref_site[1];
      delta[2] = site[2] - ref_site[2];
      residual_sum += ( (delta[0]*delta[0]*weight)+
                        (delta[1]*delta[1]*weight)+
                        (delta[2]*delta[2]*weight) );
      gradient[0] = delta[0]*2.0*weight;
      gradient[1] = delta[1]*2.0*weight;
      gradient[2] = delta[2]*2.0*weight;
      gradient_array[ i_seqs[0] ] += gradient;
    }
    //std::cout << "reference coord residual: " << residual_sum << std::endl;
    return residual_sum;
  }
}}

 /*target += (d[0]**2*w+d[1]**2*w+d[2]**2*w)
        if(gradient_array is not None):
          a = (d[0]*2*w, d[1]*2*w, d[2]*2*w)
          b = gradient_array[i_seq]
          r = ((a[0]+b[0]), (a[1]+b[1]), (a[2]+b[2]))
          gradient_array[i_seq] = r
        cntr += 1*/

