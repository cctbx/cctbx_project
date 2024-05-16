#ifndef SCITBX_MATH_G_FUNCTION_H
#define SCITBX_MATH_G_FUNCTION_H

#include <scitbx/array_family/shared.h>
#include <scitbx/constants.h>

namespace scitbx { namespace math {

//! G-function: exact and fast using look-up table.

namespace g_function {

template <typename FloatType>
FloatType Gfunction(FloatType twoPiRS)
{
  // G-function (Fourier transform of a sphere of radius r at resolution s).
  static FloatType EPS(0.001);
  if(twoPiRS > EPS) {
    return 3*(std::sin(twoPiRS)-twoPiRS*std::cos(twoPiRS))/fn::pow3(twoPiRS);
  }
  else {
    if(twoPiRS > 0) return 1-fn::pow2(twoPiRS)/10;
    else            return 1;
  }
}

template <typename FloatType>
af::shared<FloatType> Gfunction(af::const_ref<FloatType> const& s,
                                FloatType const& R,
                                bool const& volume_scale)
{
  // G-function (Fourier transform of a sphere of radius R at resolution s).
  static FloatType EPS(0.001);
  af::shared<FloatType> r(s.size());
  FloatType volume=1;
  if(volume_scale) volume = 4*scitbx::constants::pi*fn::pow3(R)/3;
  for(std::size_t i = 0; i < s.size(); i++) {
    FloatType twoPiRS = scitbx::constants::two_pi*s[i]*R;
    r[i] = volume * Gfunction(twoPiRS);
  }
  return r;
}

template <typename FloatType>
FloatType GfuncOfRSsqr(FloatType rsSqr)
{
  if(rsSqr<0) {
    if(std::abs(rsSqr)<1.e-9) rsSqr=0;
    SCITBX_ASSERT(rsSqr>=0);
  }
  return Gfunction(scitbx::constants::two_pi*std::sqrt(rsSqr));
}

template <typename FloatType>
FloatType dGfunc_by_dR(FloatType r, FloatType s)
{
  static FloatType EPS(0.001);
  FloatType twoPiRS(scitbx::constants::two_pi*r*s);
  if (std::abs(twoPiRS) > EPS)
    return 3*(3*twoPiRS*std::cos(twoPiRS) +
      (scitbx::fn::pow2(twoPiRS)-3.)*std::sin(twoPiRS)) /
      (r*scitbx::fn::pow3(twoPiRS));
  else
    return -scitbx::fn::pow2(twoPiRS)/(5*r);
}

template <typename FloatType>
af::shared<std::pair<FloatType,FloatType> >
getGfuncOfRSsqrTable(int tbllen, FloatType argmax)
{
  af::shared<std::pair<FloatType,FloatType> > GfuncOfRSsqrTable;
  GfuncOfRSsqrTable.resize(tbllen+1);
  FloatType divarg(argmax/tbllen);
  for (unsigned i = 0; i < tbllen+1; i++)
  {
    FloatType x0(GfuncOfRSsqr((i  )*divarg));
    FloatType x1(GfuncOfRSsqr((i+1)*divarg));
    GfuncOfRSsqrTable[i].first = x0;
    GfuncOfRSsqrTable[i].second = x1 - x0;
  }
  return GfuncOfRSsqrTable;
}

template <typename FloatType>
FloatType GfuncOfRSsqr_approx(FloatType arg)
{
  // Accuracy about 1 part in 10^8 or better, zero to first root (return zero
  // for larger args)
  // Quadratic interpolation would use smaller tables, but is slower
  static int tbllen(8192);
  static FloatType argcut(0.51143712);
  static FloatType argfac(tbllen/argcut);
  static af::shared<std::pair<FloatType,FloatType> > GfuncOfRSsqrTable =
    getGfuncOfRSsqrTable(tbllen,argcut);
  FloatType fracTableArg;
  int index;
  if (arg < argcut)
  {
    fracTableArg = argfac*arg;
    index = static_cast<int>(fracTableArg);
    return GfuncOfRSsqrTable[index].first +
      (fracTableArg-index)*GfuncOfRSsqrTable[index].second;
  }
  else
    return 0;
}


}}} // namespace scitbx::math::g_function

#endif // SCITBX_MATH_G_FUNCTION_H
