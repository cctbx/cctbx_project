#include <scitbx/constants.h>
#include <boost/math/special_functions/bessel.hpp>

/*
  Adopted from Randy Read's code by Pavel Afonine, 12-MAR-2014
*/

namespace scitbx { namespace math {
namespace parabolic_cylinder_d {

double dvsa(double,double);
double dvla(double,double);
double vvla(double,double);

double dv(double va, double x)
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
  if (ax <= 5.8)
    pd = dvsa(va,x);
  else
    pd = dvla(va,x);
  return pd;
}

double dvsa(double va, double x)
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

double dvla(double va,double x)
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

double vvla(double va,double x)
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
