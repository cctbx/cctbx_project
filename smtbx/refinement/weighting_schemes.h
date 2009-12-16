/// Weighting schemes for L.S. (least_squares.h)

#include <scitbx/array_family/shared.h>
#include <smtbx/import_scitbx_af.h>

#include <algorithm>
#include <cmath>

namespace smtbx { namespace refinement { namespace least_squares {

  template<typename T, class WeightingScheme>
  af::shared<T> weights(WeightingScheme const &weighting_scheme,
                        af::const_ref<T> const &fo_sq,
                        af::const_ref<T> const &sigmas,
                        af::const_ref<T> const &fc_sq)
  {
    af::shared<T> result(fo_sq.size());
    for (int i=0; i<fo_sq.size(); i++) {
      result[i] = weighting_scheme(fo_sq[i], sigmas[i], fc_sq[i]);
    }
    return result;
  }

  template <typename T>
  struct unit_weighting
  {
    T operator()(T fo_sq, T sigma, T fc_sq) const {
      return 1.;
    }
  };


  template <typename T>
  struct mainstream_shelx_weighting
  {
    T a, b;

    mainstream_shelx_weighting(T a=0.1, T b=0)
      : a(a), b(b)
    {}

    T operator()(T fo_sq, T sigma, T fc_sq) const {
      T p = (std::max(fo_sq, 0.) + 2*fc_sq)/3.;
      return 1./(sigma*sigma + std::pow(a*p, 2) + b*p);
    }
  };


}}}
