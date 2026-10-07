#ifndef SCITBX_MATH_PARABOLIC_CYLINDER_D_H
#define SCITBX_MATH_PARABOLIC_CYLINDER_D_H

#include <scitbx/constants.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/erf.hpp>

/*
  Adopted from Randy Read's code by Pavel Afonine, 12-MAR-2014
*/

namespace scitbx { namespace math {
namespace parabolic_cylinder_d {

inline double dvsa(double,double);
inline double dvla(double,double);
inline double vvla(double,double);

//! Dv(x) from modified Bessel functions or erfc, for va = -1/2, -1, -3/2
//! (x > 0) and va = -1/2 (x < 0); returns false for other cases.
/*! Used where the series (dvsa) and asymptotic (dvla) expansions both lose
    accuracy, near their crossover at |x| = 5.8 (relative error up to 2e-7
    for va = -3/2). With z = x^2/4:
      D_{-1/2}(x) = sqrt(x/(2 pi)) K_{1/4}(z)
      D_{1/2}(x)  = x^{3/2}/(2 sqrt(2 pi)) (K_{1/4}(z) + K_{3/4}(z))
      D_{-3/2}(x) = 2 (D_{1/2}(x) - x D_{-1/2}(x))   (recurrence)
      D_{-1}(x)   = sqrt(pi/2) exp(z) erfc(x/sqrt(2))
      D_{-1/2}(-x) = sqrt(pi x)/2 (I_{-1/4}(z) + I_{1/4}(z))
 */
inline bool dv_closed_form(double va, double x, double& pd)
{
  const double pi = scitbx::constants::pi;
  const double z = x * x / 4.;
  if (x > 0.) {
    if (va == -0.5 || va == -1.5) {
      const double k14 = boost::math::cyl_bessel_k(0.25, z);
      const double dm05 = std::sqrt(x / (2. * pi)) * k14;
      if (va == -0.5) { pd = dm05; return true; }
      const double dp05 = std::pow(x, 1.5) / (2. * std::sqrt(2. * pi))
                        * (k14 + boost::math::cyl_bessel_k(0.75, z));
      pd = 2. * (dp05 - x * dm05);
      return true;
    }
    if (va == -1.) {
      pd = std::sqrt(pi / 2.) * std::exp(z)
         * boost::math::erfc(x / std::sqrt(2.));
      return true;
    }
  }
  else if (x < 0. && va == -0.5) {
    const double y = -x;
    pd = std::sqrt(pi * y) / 2.
       * (boost::math::cyl_bessel_i(-0.25, z) + boost::math::cyl_bessel_i(0.25, z));
    return true;
  }
  return false;
}

inline double dv(double va, double x)
{
/* Compute parabolic cylinder function Dv(x)
   Equivalent to Mathematica ParabolicCylinderD[va,x]

   Derived from routines in program mpbdv.for by Shanjie Zhang and Jianming Jin
   distributed with their book "Computation of Special Functions",
   Copyright 1996 by John Wiley & Sons, Inc.
   Permission has been granted for purchasers of the book to incorporate these
   programs as long as the copyright is acknowledged.

   va = order of the parabolic cylinder function Dv
   x  = argument
*/
  double ax(std::abs(x));
  double pd;
  if (ax > 3.5 && ax < 9. && dv_closed_form(va, x, pd))
    return pd;
  if (ax <= 5.8)
    pd = dvsa(va,x);
  else
    pd = dvla(va,x);
  return pd;
}

inline double dvsa(double va, double x)
{
// Compute parabolic cylinder function Dv(x) for small values of |x| (<=5.8)
  static double EPS(std::pow(10.,-15));
  static double SQRT2(std::sqrt(2.));
  static double SQRTPI(std::sqrt(scitbx::constants::pi));
  double pd;
  double ep = exp(-0.25*x*x);
  double va0 = 0.5*(1.-va);
  if (va == 0.)
    pd = ep;
  else
  {
    if (x == 0.)
    {
      double ftol = std::numeric_limits<double>::epsilon();
      if (va0 <= 0. && std::abs(va0-std::floor(va0+0.5)) < ftol) pd = 0.;
      else pd = SQRTPI/(boost::math::tgamma(va0)*std::pow(2,-0.5*va));
    }
    else
    {
      double a0 = std::pow(2,-0.5*va-1.)*ep/boost::math::tgamma(-va);
      double vt = -0.5*va;
      pd = boost::math::tgamma(vt);
      double r(1.),r1(pd);
      int m(1);
      while (m<=250 && std::abs(r1)>=std::abs(pd)*EPS)
      {
        double vm = 0.5*(m-va);
        r = -r*SQRT2*x/m;
        r1 = boost::math::tgamma(vm)*r;
        pd += r1;
        m++;
      }
      pd *= a0;
    }
  }
  return pd;
}

inline double dvla(double va,double x)
{
// Compute parabolic cylinder function Dv(x) for large values of |x| (>5.8)
  static double EPS(std::pow(10.,-12));
  double xsqr(x*x);
  double ep(exp(-0.25*xsqr));
  double a0(std::pow(std::abs(x),va)*ep);
  double r(1.);
  double pd(1.);
  int k(1);
  while (k<=16 && std::abs(r/pd)>=EPS)
  {
    r = -0.5*r*(2.*k-va-1.)*(2.*k-va-2.)/(k*xsqr);
    pd += r;
    k++;
  }
  pd *= a0;
  if (x < 0.)
  {
    double x1(-x);
    double vl = vvla(va,x1);
    pd = scitbx::constants::pi*vl/boost::math::tgamma(-va) +
         cos(scitbx::constants::pi*va)*pd;
  }
  return pd;
}

inline double vvla(double va,double x)
{
// Compute parabolic cylinder function Vv(x) for large argument
  static double EPS(std::pow(10.,-12));
  static double SQRT2BYPI(std::sqrt(2./scitbx::constants::pi));
  double xsqr(x*x);
  double qe(exp(0.25*xsqr));
  double a0 = std::pow(std::abs(x),-va-1.)*SQRT2BYPI*qe;
  double r(1.),pv(1.);
  int k(1);
  while (k<=18 && std::abs(r/pv)>=EPS)
  {
    r = 0.5*r*(2.*k+va-1.)*(2.*k+va)/(k*xsqr);
    pv += r;
    k++;
  }
  pv *= a0;
  if (x < 0.)
  {
    double x1(-x);
    double pdl = dvla(va,x1);
    double gl = boost::math::tgamma(-va);
    double piva = scitbx::constants::pi*va;
    double dsl = std::pow(sin(piva),2);
    pv = dsl*gl/scitbx::constants::pi*pdl - cos(piva)*pv;
  }
  return pv;
}

}}}

#endif // SCITBX_MATH_PARABOLIC_CYLINDER_D_H
