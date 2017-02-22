#ifndef MMTBX_RIDING_H
#define MMTBX_RIDING_H

#include <iostream>
#include <string>
#include <cctbx/sgtbx/space_group.h>
#include <boost/python/list.hpp>
//#include <assert.h>
#include <mmtbx/error.h>

using namespace std;
namespace mmtbx { namespace hydrogens {
namespace af=scitbx::af;
using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;


//template <typename FloatType=double>
class riding_coefficients
{
  public:
    string htype;
    int ih;
    int a0;
    int a1;
    int a2;
    int a3;
    double a;
    double b;
    double h;
    int n;
    double disth;
    riding_coefficients() {
    }
    riding_coefficients(
         string const& htype_,
         int const& ih_,
         int const& a0_,
         int const& a1_,
         int const& a2_,
         int const& a3_,
         double const& a_,
         double const& b_,
         double const& h_,
         int const& n_,
         double const& disth_)
    :
      htype(htype_), ih(ih_), a0(a0_), a1(a1_), a2(a2_), a3(a3_), a(a_), b(b_),
      h(h_), n(n_), disth(disth_)
    {}
};

vec3<double> compute_h_position(
                   riding_coefficients rc,
                   af::shared<vec3<double> > const& sites_cart)
  {
    vec3<double> r0 = sites_cart[rc.a0];
    vec3<double> rh = sites_cart[rc.ih];
    vec3<double> r1 = sites_cart[rc.a1];
    vec3<double> rh_calc;
    double a = rc.a;
    double b = rc.b;
    double h = rc.h;
    double dh = rc.disth;
    if (rc.htype == "flat_2neigbs") {
      vec3<double> r2 = sites_cart[rc.a2];
      vec3<double> u10 = (r1-r0).normalize();
      vec3<double> u20 = (r2-r0).normalize();
      vec3<double> rh0 = a * u10 + b * u20;
      double rh0_length = rh0.length();
      MMTBX_ASSERT(rh0_length > 0.);
      vec3<double> uh0 = rh0 / rh0_length;
      rh_calc = r0 + dh * uh0;
    } else if (rc.htype == "2neigbs") {
      vec3<double> r2 = sites_cart[rc.a2];
      vec3<double> u10 = (r1-r0).normalize();
      vec3<double> u20 = (r2-r0).normalize();
      vec3<double> v0 = (u10.cross(u20)).normalize();
      vec3<double> rh0 = (a * u10 + b * u20 + h * v0);
      double rh0_length = rh0.length();
      MMTBX_ASSERT(rh0_length > 0.);
      vec3<double> uh0 = rh0 / rh0_length;
      rh_calc = r0 + dh * uh0;
    } else if (rc.htype == "2tetra") {
      vec3<double> r2 = sites_cart[rc.a2];
      vec3<double> u10 = (r1-r0).normalize();
      vec3<double> u20 = (r2-r0).normalize();
      vec3<double> v0 = (u10.cross(u20)).normalize();
      vec3<double> d0 = (a * u10 + b * u20).normalize();
      rh_calc = r0 + dh * (cos(h) * d0 + sin(h) * v0);
    } else if (rc.htype == "3neigbs") {
      vec3<double> r2 = sites_cart[rc.a2];
      vec3<double> r3 = sites_cart[rc.a3];
      vec3<double> u10 = (r1-r0).normalize();
      vec3<double> u20 = (r2-r0).normalize();
      vec3<double> u30 = (r3-r0).normalize();
      vec3<double> rh0 = (a * u10 + b * u20 + h * u30);
      double rh0_length = rh0.length();
      MMTBX_ASSERT(rh0_length > 0.);
      vec3<double> uh0 = rh0 / rh0_length;
      rh_calc = r0 + dh * uh0;
    } else if (rc.htype=="alg1b" || rc.htype=="alg1a" || rc.htype=="prop") {
      //std::cout << "Type 1neigb\n";
      vec3<double> rb1 = sites_cart[rc.a2];
      int n = rc.n;
      double phi = b + n * 2.0 * scitbx::constants::pi/3.0;
      double alpha = a;
      double salpha = sin(alpha);
      double calpha = cos(alpha);
      double sphi = sin(phi);
      double cphi = cos(phi);
      vec3<double> u1 = (r0 - r1).normalize();
      vec3<double> rb10 = rb1 -r1;
      vec3<double> u2 = (rb10 - (rb10*u1) * u1).normalize();
      vec3<double> u3 = u1.cross(u2);
      rh_calc = r0 + dh * (salpha*(cphi*u2 + sphi*u3) - calpha*u1);
    } else {
      rh_calc = rh;
    }

    return rh_calc;
} // compute_H_position

// returns new sites_cart (does not overwrite input sites_cart)
af::shared<vec3<double> > apply_new_H_positions(
    af::shared<vec3<double> > const& sites_cart,
    boost::python::list const& parameterization_)
  {
// make local list parameterization
    af::shared<riding_coefficients> parameterization;
    for(std::size_t i=0;i<boost::python::len(parameterization_);i++) {
      parameterization.push_back(
          boost::python::extract<riding_coefficients>(parameterization_[i])());
    }
// create new sites_cart and initialize
    af::shared<vec3<double> > sites_cart_new(sites_cart.size());
    for(std::size_t i=0; i < sites_cart.size(); i++) {
      sites_cart_new[i] = sites_cart[i];
    }
// calculate positions of H atoms according to parameterization
    for(int i=0; i < parameterization.size(); i++) {
      riding_coefficients rc = parameterization[i];
      int ih = rc.ih;
      vec3<double> rh_calc;
      rh_calc = compute_h_position(rc, sites_cart);
      sites_cart_new[ih] = rh_calc;
    }
  return sites_cart_new;
}

vec3<double> G_unitvector(
                   vec3<double> Gu,
                   vec3<double> r)
{
  vec3<double> Gr;
  double Gx, Gy, Gz;
  double x = r[0], y = r[1], z = r[2];
  double Gxu = Gu[0], Gyu = Gu[1], Gzu = Gu[2];
  double denom = pow(r.length(),3);
  MMTBX_ASSERT(denom > 0.);
  Gx = ((pow(y,2) + pow(z,2)) * Gu[0] - x*y*Gu[1] - x*z*Gu[2]) / denom;
  Gy = ( -x*y*Gu[0] + (pow(x,2) + pow(z,2))*Gu[1] - y*z*Gu[2]) / denom;
  Gz = ( -x*z*Gu[0] - y*z*Gu[1] + (pow(x,2) + pow(y,2))*Gu[2]) / denom;
  Gr[0] = Gx;
  Gr[1] = Gy;
  Gr[2] = Gz;
  return Gr;
}

// transformation for the derivative of a cross product
af::shared<vec3<double> > G_crossproduct(
                   vec3<double> Gv,
                   vec3<double> r1,
                   vec3<double> r2)
{
  af::shared<vec3<double> > G_r1_r2;
  vec3<double> Gr1, Gr2;
  double x1 = r1[0], y1 = r1[1], z1 = r1[2];
  double x2 = r2[0], y2 = r2[1], z2 = r2[2];
  double Gx_v = Gv[0], Gy_v = Gv[1], Gz_v = Gv[2];
  Gr1[0] = -z2 * Gy_v + y2 * Gz_v;
  Gr1[1] =  z2 * Gx_v - x2 * Gz_v;
  Gr1[2] = -y2 * Gx_v + x2 * Gy_v;
  Gr2[0] =  z1 * Gy_v - y1 * Gz_v;
  Gr2[1] = -z1 * Gx_v + x1 * Gz_v;
  Gr2[2] =  y1 * Gx_v - x1 * Gy_v;
  G_r1_r2[0] = Gr1;
  G_r1_r2[1] = Gr2;
  return G_r1_r2;
}

}} // namespace mmtbx::riding_h

#endif // MMTBX_RIDING_H
