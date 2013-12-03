#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <cmath>

namespace scitbx { namespace math {
namespace gaussian_fit_1d_analytical {

//! Analytical one-Gaussian a*exp(-b*x**2) approximation of y = f(x).
template <typename FloatType=double>
class compute
{
public:
  FloatType a, b;

  compute() {}

  compute(
    af::const_ref<FloatType> const& x,
    af::const_ref<FloatType> const& y)
  :
  a(0), b(0)
  {
    SCITBX_ASSERT(x.size() == y.size());
    FloatType p=0., q=0., r=0., s=0.;
    for(std::size_t i=0; i < y.size(); i++) {
      FloatType d_ = y[i];
      if(d_<=0) goto done;
      FloatType d = std::log(d_);
      FloatType v = x[i]*x[i];
      p += d;
      q += v;
      r += (v*v);
      s += (v*d);
    }
    if(r != 0) {
      int n = y.size();
      FloatType den = n-q*q/r;
      if(den != 0) {
        FloatType u = (p-s*q/r)/den;
        b = (u*q-s)/r;
        a = std::exp(u);
      }
    }
  done: ;
  }

  compute(
    af::const_ref<FloatType> const& x,
    af::const_ref<FloatType> const& y,
    af::const_ref<FloatType> const& z)
  :
  a(0), b(0)
  {
    SCITBX_ASSERT(x.size() == y.size());
    SCITBX_ASSERT(x.size() == z.size());
    FloatType p=0., q=0., r=0., s=0.;
    for(std::size_t i=0; i < y.size(); i++) {
      FloatType zi_ = z[i];
      if(zi_==0) goto done;
      FloatType d_ = y[i]/zi_;
      if(d_<=0) goto done;
      FloatType d = std::log(d_);
      FloatType v = x[i]*x[i];
      p += d;
      q += v;
      r += (v*v);
      s += (v*d);
    }
    if(r != 0) {
      int n = y.size();
      FloatType den = n-q*q/r;
      if(den != 0) {
        FloatType u = (p-s*q/r)/den;
        b = (u*q-s)/r;
        a = std::exp(u);
      }
    }
  done: ;
  }
};
}}} // namespace scitbx::math::gaussian_fit_1d_analytical
