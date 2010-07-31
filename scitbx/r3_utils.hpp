#ifndef SCITBX_R3_UTILS_HPP
#define SCITBX_R3_UTILS_HPP

#include <scitbx/vec3.h>
#include <tbxx/error_utils.hpp>
#include <boost/unordered_set.hpp>
#include <boost/noncopyable.hpp>
#include <vector>

namespace scitbx { namespace r3_utils {

  struct clash_detector_simple : boost::noncopyable
  {
    std::vector<boost::unordered_set<unsigned> > exclusion_sets;
    double threshold_sq;

    clash_detector_simple() {}

    clash_detector_simple(
      unsigned n_vertices,
      double threshold)
    :
      exclusion_sets(n_vertices),
      threshold_sq(threshold*threshold)
    {}

    void
    add_exclusion(
      int i,
      int j)
    {
      TBXX_ASSERT(i < j);
      exclusion_sets[i].insert(j);
    }

    bool
    has_clash(
      af::const_ref<vec3<double> > const& sites_cart)
    {
      TBXX_ASSERT(sites_cart.size() == exclusion_sets.size());
      unsigned n = static_cast<unsigned>(sites_cart.size());
      for(unsigned i=0;i<n-1;i++) {
        vec3<double> const& sci = sites_cart[i];
        boost::unordered_set<unsigned> const& esi = exclusion_sets[i];
        boost::unordered_set<unsigned>::const_iterator
          esi_end = esi.end();
        for(unsigned j=i+1;j<n;j++) {
          if (esi.find(j) != esi_end) continue;
          vec3<double> const& scj = sites_cart[j];
          double d = (scj - sci).length_sq();
          if (d < threshold_sq) {
            return true;
          }
        }
      }
      return false;
    }
  };

}}

#endif // GUARD
