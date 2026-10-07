#ifndef CCTBX_FRENCH_WILSON_H
#define CCTBX_FRENCH_WILSON_H

#include <algorithm>
#include <cmath>
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

//! 2nd/4th-moment Rice approximation to the French-Wilson posterior of E.
/*! Effective Rice parameters (Dobs, Eeff) for an observed intensity with
    measurement error, matching the 2nd and 4th moments of the Rice
    distribution for E to those of the French-Wilson posterior.

    Writing J = E^2, the posterior for J is proportional to
      acentric: exp(-(J-m)^2/(2 sigesq^2))           m = eosq - sigesq^2
      centric:  J^(-1/2) exp(-(J-m)^2/(2 sigesq^2))  m = eosq - sigesq^2/2
    on J >= 0. Integrating J^k (J-m) times this by parts gives
      <J^2> = m <J> + k sigesq^2,   k = 1 (acentric), 1/2 (centric)
    so <E^4> follows exactly from the French-Wilson <E^2>.

    The Rice moments are <E^2> = 1 + eta and
      acentric: <E^4> = 2 + 4 eta + 2 eta^2 -   (D^2 + eta)^2
      centric:  <E^4> = 3 + 6 eta + 3 eta^2 - 2 (D^2 + eta)^2
    with eta = D^2 (Eeff^2 - 1). The matching solution is D^2 = S - eta with
      acentric: S = sqrt(2 eta^2 + 4 eta - zeta),      zeta = <E^4> - 2
      centric:  S = sqrt((3 eta^2 + 6 eta - zeta)/2),  zeta = <E^4> - 3
    and Eeff^2 = 1 + eta/D^2. There is no real solution (valid = false) if
    the discriminant is negative or D^2 <= 0: the posterior is more
    dispersed than any Rice distribution with the same <E^2> (in practice
    very strong but very imprecise intensities); approaching that limit,
    D -> 0 and Eeff grows without bound. sigesq <= 0 is taken as error-free
    data (D = 1).
 */
struct rice_moments
{
  bool   valid;
  double mu2;   // posterior <E^2>
  double mu4;   // posterior <E^4>
  double dsqr;  // Dobs^2
  double dobs;
  double eeff;
  rice_moments()
  : valid(false), mu2(0), mu4(0), dsqr(0), dobs(0), eeff(0) {}
};

inline
rice_moments
rice_from_intensity(double eosq, double sigesq, bool centric)
{
  rice_moments result;
  if (!(sigesq > 0.)) {
    result.valid = true;
    result.mu2 = std::max(eosq, 0.);
    result.mu4 = result.mu2 * result.mu2;
    result.dsqr = result.dobs = 1.;
    result.eeff = std::sqrt(result.mu2);
    return result;
  }
  const double k = centric ? 0.5 : 1.;
  const double varobs = sigesq * sigesq;
  const double m = eosq - k * varobs;
  // the closed form avoids the cancellation in m + sigesq*lambda for m << 0,
  // which matters near the validity boundary where D^2 is ill-conditioned
  result.mu2 = expectEsqFW(eosq, sigesq, centric);
  result.mu4 = m * result.mu2 + k * varobs;
  const double eta = result.mu2 - 1.;
  double gap, disc;
  if (centric) {
    const double zeta = result.mu4 - 3.;
    gap = eta * eta + 6. * eta - zeta;   // 2(S^2 - eta^2)
    disc = (3. * eta * eta + 6. * eta - zeta) / 2.;
  }
  else {
    const double zeta = result.mu4 - 2.;
    gap = eta * eta + 4. * eta - zeta;   // S^2 - eta^2
    disc = 2. * eta * eta + 4. * eta - zeta;
  }
  if (!(disc >= 0.)) return result;      // also catches NaN
  const double S = std::sqrt(disc);
  // S - eta suffers cancellation when eta > 0 and S ~ eta: rationalize
  double dsqr = (eta >= 0.) ? (centric ? gap / 2. : gap) / (eta + S)
                            : S - eta;
  if (!(dsqr > 0.)) return result;
  dsqr = std::min(dsqr, 1.);
  result.valid = true;
  result.dsqr = dsqr;
  result.dobs = std::sqrt(dsqr);
  result.eeff = std::sqrt(std::max(1. + eta / dsqr, 0.));
  return result;
}

//! Dobs given to reflections with no Rice solution: negligible weight in
//! the Rice approximation (they need the exact likelihood).
inline double rice_dobs_none() { return 0.01; }

//! (Dobs, Eeff) for a reflection: the moment-matched values, or
//! (rice_dobs_none(), sqrt(<E^2>)) where there is no Rice solution, so that
//! Feff stays a sensible amplitude.
inline
void
rice_dobs_eeff(rice_moments const& rice, double& dobs, double& eeff)
{
  if (rice.valid) {
    dobs = rice.dobs;
    eeff = rice.eeff;
  }
  else {
    dobs = rice_dobs_none();
    eeff = std::sqrt(std::max(rice.mu2, 0.));
  }
}

//! Mean and standard deviation of sqrt(u) for the French-Wilson posterior
//! u^p exp(-(u - h)^2/2) on u >= 0 (p = 0 acentric, -1/2 centric): the
//! closed forms with sigesq = 1 and eosq = h + 1 (acentric) or h + 1/2
//! (centric).
inline
void
french_wilson_sqrt_u_moments(double h, bool centric, double& mean, double& sd)
{
  const double eosq = h + (centric ? 0.5 : 1.0);
  mean = expectEFW(eosq, 1.0, centric);
  const double mean_u = expectEsqFW(eosq, 1.0, centric);
  sd = std::sqrt(std::max(mean_u - mean * mean, 0.0));
}

//! Recover (I, SIGI) from French-Wilson (F, SIGF), given the prior mean
//! intensity the French-Wilson calculation used.
/*! In units of SIGI the posterior of u = J/SIGI is u^p exp(-(u - h)^2/2),
    with h = I/SIGI - SIGI/<I> (acentric) or I/SIGI - SIGI/(2<I>)
    (centric), so F = <sqrt(u)> sqrt(SIGI) and SIGF = sd(sqrt(u)) sqrt(SIGI).
    SIGF/F depends on h alone, decreasing monotonically from sqrt(4/pi - 1)
    (acentric) or sqrt(pi/2 - 1) (centric) as h -> -infinity to 1/(2h) for
    large h, so h follows from SIGF/F by bisection; then
    SIGI = (F/<sqrt(u)>)^2 and I = SIGI*(h + c*SIGI/<I>), c = 1 (acentric)
    or 1/2 (centric). mean_intensity only enters the last, weak-data term.
    Where SIGF/F is at or above its value at h = h_min, h is set to h_min
    and prior_dominated is true: such amplitudes carry essentially no
    information about I. Returns false (I, SIGI, h unchanged) if F, SIGF or
    mean_intensity is not positive.
 */
inline
bool
invert_french_wilson(double F, double SIGF, double mean_intensity,
  bool centric, double& I, double& SIGI, double& h, bool& prior_dominated,
  double h_min = -6.0)
{
  if (!(F > 0 && SIGF > 0 && mean_intensity > 0)) return false;
  const double r = SIGF / F;
  double mean, sd;
  french_wilson_sqrt_u_moments(h_min, centric, mean, sd);
  h = h_min;
  prior_dominated = !(r < sd / mean);
  if (!prior_dominated) {
    double lo = h_min, hi = std::max(10.0, 2.0 / r);
    for (int it = 0; it < 60; it++) {
      french_wilson_sqrt_u_moments(hi, centric, mean, sd);
      if (sd / mean < r) break;
      hi *= 2;
    }
    for (int it = 0; it < 80; it++) {
      double mid = 0.5 * (lo + hi);
      french_wilson_sqrt_u_moments(mid, centric, mean, sd);
      if (sd / mean > r) lo = mid; else hi = mid;
      if (hi - lo < 1e-10 * std::max(1.0, std::abs(hi))) break;
    }
    h = 0.5 * (lo + hi);
  }
  french_wilson_sqrt_u_moments(h, centric, mean, sd);
  const double sqrt_sigi = F / mean;
  SIGI = sqrt_sigi * sqrt_sigi;
  const double c = centric ? 0.5 : 1.0;
  I = SIGI * (h + c * SIGI / mean_intensity);
  return true;
}

} // namespace cctbx

#endif // CCTBX_FRENCH_WILSON_H
