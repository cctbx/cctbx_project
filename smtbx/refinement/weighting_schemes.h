/// Weighting schemes for L.S. (least_squares.h)

#include <algorithm>
#include <cmath>

namespace smtbx { namespace refinement { namespace least_squares {

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
      return 1./(sigma + std::pow(a*p, 2) + b*p);
    }
  };


}}}
