#include <scitbx/constants.h>
#include <boost/math/special_functions/bessel.hpp>
#include <scitbx/math/parabolic_cylinder_d.h>
#include <scitbx/math/erf.h>
#include <cctbx/error.h>

namespace cctbx {
namespace pc=scitbx::math::parabolic_cylinder_d;
namespace fn=scitbx::fn;

template <typename floatType>
floatType expectEFWacen(floatType eosq, floatType sigesq)
{
/* Acentric: Compute French & Wilson posterior expected value of E, from the
   normalised observed intensity (Eobs^2=Iobs/<I>) and its standard deviation
*/
  const floatType CROSSOVER1(-12.5), CROSSOVER2(18.);
  static floatType SQRT2(std::sqrt(2.));
  floatType ee;
  floatType x((eosq-fn::pow2(sigesq))/sigesq);
  floatType xsqr(fn::pow2(x));
  if (x < CROSSOVER1) // Large negative argument: asymptotic approximation
    ee = std::sqrt(-scitbx::constants::pi*sigesq/x) *
         (-916620705. + xsqr *
         (91891800.   + xsqr *
         (-11531520.  + xsqr *
         (1935360.    + xsqr *
         (-491520.    + xsqr * 262144.))))) /
         (-495452160. + xsqr *
         (55050240.   + xsqr *
         (-7864320.   + xsqr *
         (1572864.    + xsqr *
         (-524288.    + xsqr * 524288.)))));
  else if (x > CROSSOVER2) // Large positive argument: asymptotic approximation
    ee = std::sqrt(sigesq) *
         (-45045. + 32.*xsqr *
         (-315.    + 8.*xsqr *
         (-15.    - 16.*xsqr + 128.*fn::pow2(xsqr)))) /
         (32768*std::pow(x,7.5));
  else // Moderate arguments: analytical integral
    ee = std::sqrt(sigesq/2.) * std::exp(-xsqr/4.) *
         pc::dv(-1.5,-x) / scitbx::math::erfc(-x/SQRT2);
  return ee;
}

template <typename floatType>
floatType expectEsqFWacen(floatType eosq, floatType sigesq)
{
/* Acentric: Compute French & Wilson posterior expected value of E^2, from the
   normalised observed intensity (Eobs^2=Iobs/<I>) and its standard deviation
*/
  const floatType CROSSOVER1(-8.9), CROSSOVER2(5.7);
  static floatType SQRT2BYPI(std::sqrt(2./scitbx::constants::pi));
  static floatType SQRT2(std::sqrt(2.));
  floatType eesq((eosq-fn::pow2(sigesq))); // Default for significantly positive x
  floatType x(eesq/(SQRT2*sigesq));
  floatType xsqr(fn::pow2(x));
  if (x < CROSSOVER1) // Large negative argument: asymptotic approximation
    eesq *= (-135135 + xsqr *
            (20790   + xsqr *
            (-3780   + xsqr *
            (840     + xsqr *
            (-240    + xsqr *
            (96      - xsqr * 64)))))) /
            (-135135 + xsqr *
            (20790   + xsqr *
            (-3780   + xsqr *
            (840     + xsqr *
            (-240    + xsqr *
            (96      + xsqr *
            (-64     + xsqr * 128)))))));
  else if (x <= CROSSOVER2) // Moderate arguments: analytical integral
    eesq += SQRT2BYPI * sigesq / (std::exp(xsqr) * scitbx::math::erfc(-x));
  return eesq;
}

template <typename floatType>
floatType expectEFWcen(floatType eosq, floatType sigesq)
{
/* Centric: Compute French & Wilson posterior expected value of E, from the
   normalised observed intensity (Eobs^2=Iobs/<I>) and its standard deviation
*/
  const floatType CROSSOVER1(-17.5), CROSSOVER2(17.5);
  static floatType SQRTPI(std::sqrt(scitbx::constants::pi));
  floatType pcdratio,ee;
  floatType x(sigesq/2.-eosq/sigesq);
  floatType xsqr(fn::pow2(x));
  if (x < CROSSOVER1) // Large negative argument: asymptotic approximation
    pcdratio = (1024.*SQRTPI*std::pow(-x,6.5)) /
               (3465. + xsqr *
               (840.  + xsqr *
               (384.  + xsqr * 1024.)));
  else if (x > CROSSOVER2) // Large positive argument: asymptotic approximation
    pcdratio = (3440640. + xsqr *
               (-491520. + xsqr *
               (98304.   + xsqr *
               (-32768.  + xsqr * 32768.)))) /
               (675675.  + xsqr *
               (-110880. + xsqr *
               (26880.   + xsqr *
               (-12288.  + xsqr * 32768.)))) / std::sqrt(x);
  else // Moderate arguments: analytical integral
    pcdratio = pc::dv(-1.,x) / pc::dv(-0.5,x);
  ee = std::sqrt(sigesq/scitbx::constants::pi)*pcdratio;
  return ee;
}

template <typename floatType>
floatType expectEsqFWcen(floatType eosq, floatType sigesq)
{
/* Centric: Compute French & Wilson posterior expected value of E^2, from the
   normalised observed intensity (Eobs^2=Iobs/<I>) and its standard deviation
*/
  const floatType CROSSOVER1(-17.5), CROSSOVER2(17.5);
  floatType pcdratio,eesq;
  floatType x(sigesq/2.-eosq/sigesq);
  floatType xsqr(fn::pow2(x));
  if (x < CROSSOVER1) // Large negative argument: asymptotic approximation
    pcdratio = (45045. + xsqr *
               (10080. + xsqr *
               (3840.  + xsqr *
               (4096.  - xsqr * 32768.)))) /
               (x *
               (55440. + xsqr *
               (13440. + xsqr *
               (6144.  + xsqr * 16384.))));
  else if (x > CROSSOVER2) // Large positive argument: asymptotic approximation
    pcdratio = (11486475. + xsqr *
               (-1441440. + xsqr *
               (241920.   + xsqr *
               (-61440.   + xsqr * 32768.)))) /
               (x *
               (675675.   + xsqr *
               (-110880.  + xsqr *
               (26880.    + xsqr *
               (-12288.   + xsqr * 32768.)))));
  else // Moderate arguments: analytical integral
    pcdratio = pc::dv(-1.5,x) / pc::dv(-0.5,x);
  eesq = sigesq*pcdratio/2.;
  return eesq;
}

template <typename floatType>
floatType expectEFW(floatType eosq, floatType sigesq, bool centric)
{
  if (sigesq <= 0.) // Apparently no measurement error
  {
    CCTBX_ASSERT(eosq>=0.); // Can only allow zero sigma for non-negative I
    return std::sqrt(eosq);
  }
  floatType eEFW;
  if (centric) eEFW = expectEFWcen(eosq,sigesq);
  else         eEFW = expectEFWacen(eosq,sigesq);
  return eEFW;
}

template <typename floatType>
floatType expectEsqFW(floatType eosq, floatType sigesq, bool centric)
{
  if (sigesq <= 0.)
  {
    CCTBX_ASSERT(eosq>=0.);
    return eosq;
  }
  floatType eEsqFW;
  if (centric) eEsqFW = expectEsqFWcen(eosq,sigesq);
  else         eEsqFW = expectEsqFWacen(eosq,sigesq);
  return eEsqFW;
}

template <typename floatType, typename af_float, typename bool1D>
bool is_FrenchWilson(af_float F, af_float SIGF, bool1D is_centric,
                     floatType eps=0.001)
{
  int nviolations(0);
  int NREFL(F.size());
  for (unsigned r = 0; r < NREFL; r++)
  {
    if (F[r] <= 0. || SIGF[r] <= 0.) return false; // French-Wilson always positive
    floatType SIGFoverF(SIGF[r]/F[r]);
    // After French-Wilson, SIGF/F has a maximum for centrics and acentrics
    if (SIGFoverF > 1) return false; // Don't expect any big violations.
    floatType maxrat = (is_centric[r]) ? 0.756 : 0.523;
    if (SIGFoverF > maxrat) nviolations ++;
  }
  floatType fracviol(floatType(nviolations)/NREFL);
  // std::cout << "Fraction of FW violations: " << fracviol << std::endl;
  if (fracviol > eps) return false;
  else return true;
}

} // namespace cctbx
