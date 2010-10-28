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
                        af::const_ref<T> const &fc_sq,
                        T scale_factor)
  {
    af::shared<T> result(fo_sq.size());
    for (int i=0; i<fo_sq.size(); i++) {
      result[i] = weighting_scheme(fo_sq[i], sigmas[i], fc_sq[i],
                                   scale_factor);
    }
    return result;
  }

  template <typename T>
  struct unit_weighting
  {
    T operator()(T fo_sq, T sigma, T fc_sq, T scale_factor) const {
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

    /* Shelxl uses fo_sq/k, sigma/k, fc_sq in its calculation of the weights
       whereas we use fo_sq, sigma and k*fc_sq. To keep our weights compatible
       with shelxl values of a and b, we must include the factor of k to the
       term (b*P) in order to make the weights linear wrt the scale factor:

         w = 1/[sigma^2 * (a*P)^2 + k*a*P],

         where P = [max(fo_sq, 0) + 2*k*fc_sq]
     */
    T operator()(T fo_sq, T sigma, T fc_sq, T scale_factor) const {
      T p = (std::max(fo_sq, 0.) + 2*scale_factor*fc_sq)/3.;
      return 1./(sigma*sigma + std::pow(a*p, 2) + b*scale_factor*p);
    }
  };


}}}
