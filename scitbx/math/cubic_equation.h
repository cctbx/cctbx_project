#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <tbxx/error_utils.hpp>
#include <boost/optional.hpp>
#include <cmath>

namespace scitbx { namespace math {
namespace cubic_equation {

//! Analytical solution of ax**3 + bx**2 + cx + d = 0.
// Returns zero in case of imaginary roots (so the name 'real')

template <typename FTW=long double,
          typename FTO=double>
class real
{
public:
  FTW A,B,D, p,a_,b_,c_,d_;
  vec3<boost::optional<FTO> > x;
  static const FTW one_over_three;

  real() {}

  real(
    FTW const& a,
    FTW const& b,
    FTW const& c,
    FTW const& d)
  :
  A(0), B(0), D(0)
  {
    FTW const eps = 1.e-9;
    SCITBX_ASSERT(a != 0.);
    a_ = static_cast<FTW>(a);
    b_ = static_cast<FTW>(b);
    c_ = static_cast<FTW>(c);
    d_ = static_cast<FTW>(d);
    p = b_/a_;
    FTW q = c_/a_;
    FTW r = d_/a_;
    FTW pp = p*p;
    A = (3.*q-pp)/3.;
    B = (2.*p*pp-9.*p*q+27.*r)/27.;
    D = (A*A*A)/27.+(B*B)/4.;
    bool flag=std::abs(A)<eps && std::abs(B)<eps && std::abs(D)<eps;
    if(D<std::abs(eps) && D<0.)       case_3();
    else if(flag)                     case_0();
    else if(D<std::abs(eps) && D>=0.) case_1();
    else if(D>0.)                     case_2();
    else throw TBXX_UNREACHABLE_ERROR();
  }

  vec3<boost::optional<FTO> > residual() {
    vec3<boost::optional<FTO> > result;
    for(std::size_t i=0; i < 3; i++) {
      if (x[i]) {
        result[i] = static_cast<FTO>(
          a_*std::pow((*x[i]),3)+b_*std::pow((*x[i]),2)+c_*(*x[i])+d_);
      }
    }
    return result;
  }

  FTW fractional_power(FTW const& arg, FTW const& pwr) {
    FTW r;
    arg<0. ? r = -std::pow(-arg, pwr) : r = std::pow(arg, pwr);
    return r;
  }

  void case_0() {
    x[0] = static_cast<FTO>(-fractional_power(d_/a_, one_over_three));
    x[1] = static_cast<FTO>(*x[0]);
    x[2] = static_cast<FTO>(*x[0]);
  }

  void case_1() {
    FTW sqrtD = std::sqrt(D);
    FTW minus_b_over_2 = -B/2.;
    FTW M = minus_b_over_2 + sqrtD;
    FTW N = minus_b_over_2 - sqrtD;
    M = fractional_power(M, one_over_three);
    N = fractional_power(N, one_over_three);
    FTW p_over_3 = p/3.;
    x[0] = static_cast<FTO>(M+N-p_over_3);
    x[1] = static_cast<FTO>(-(M+N)/2.-p_over_3);
    x[2] = *x[1];
  }

  void case_2() {
    SCITBX_ASSERT(D>=0);
    FTW sqrtD = std::sqrt(D);
    FTW minus_B_over_2 = -B/2;
    FTW R = minus_B_over_2+sqrtD;
    FTW S = fractional_power(R, one_over_three);
    FTW T = minus_B_over_2-sqrtD;
    FTW U = fractional_power(T, one_over_three);
    x[0] = static_cast<FTO>(S+U-b_/(3.*a_));
  }

  void case_3() {
    SCITBX_ASSERT(A<0.);
/* //alternative method: looks wildly different but delivers exact same results
      FTW iarg = (B*B)/4.-D;
      SCITBX_ASSERT(iarg>0.);
      FTW i = std::sqrt(iarg);
      FTW j = 0;
      if(i<0.) j = -std::pow(-i, one_over_three);
      else     j =  std::pow( i, one_over_three);
      FTW k = std::acos(-(B/(2.*i)));
      FTW k_over_3 = k/3.;
      FTW m = std::cos(k_over_3);
      FTW n = std::sqrt(3.) * std::sin(k_over_3);
      FTW p = -(b_ / (3. * a_));
      x1 = 2 * j * m + p;
      x2 = -j * (m + n) + p;
      x3 = -j * (m - n) + p;
*/
    FTW b_over_a = B/A;
    FTW b_over_a_sq = b_over_a*b_over_a;
    FTW arg = std::sqrt(-b_over_a_sq*27./(4.*A));
    if(std::abs(1-std::abs(arg))<1.e-9) arg = 1.;
    FTW phi = 0.;
    if(B>0.) phi = std::acos(-arg)/3.;
    else     phi = std::acos( arg)/3.;
    FTW ta3 = 2.*std::sqrt(-A/3.);
    FTW two_pi_over_3 = 2.*scitbx::constants::pi/3.;
    FTW p_over_3 = p/3.;
    x[0] = static_cast<FTO>(ta3*std::cos(phi              )-p_over_3);
    x[1] = static_cast<FTO>(ta3*std::cos(phi+two_pi_over_3)-p_over_3);
    x[2] = static_cast<FTO>(ta3*std::cos(phi-two_pi_over_3)-p_over_3);
  }

};

template <typename FTW, typename FTO>
const FTW real<FTW, FTO>::one_over_three = 1/3.;

}}}
