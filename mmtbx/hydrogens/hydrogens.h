#ifndef MMTBX_RIDING_H
#define MMTBX_RIDING_H

#include <iostream>
#include <string>
#include <cctbx/sgtbx/space_group.h>
#include <boost/python/list.hpp>
#include <mmtbx/error.h>

using namespace std;
namespace mmtbx { namespace hydrogens {
namespace af=scitbx::af;
using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;


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
                   riding_coefficients const& rc,
                   af::const_ref<vec3<double> > const& sites_cart)
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

//------------------------
// APPLY_NEW_H_POSITIONS
//------------------------
void apply_new_H_positions(
    af::ref<vec3<double> > sites_cart,
    boost::python::list const& parameterization_)
  {
// make local list parameterization
    af::shared<riding_coefficients> parameterization;
    for(std::size_t i=0;i<boost::python::len(parameterization_);i++) {
      parameterization.push_back(
          boost::python::extract<riding_coefficients>(parameterization_[i])());
    }
// calculate positions of H atoms according to parameterization
    for(int i=0; i < parameterization.size(); i++) {
      riding_coefficients rc = parameterization[i];
      int ih = rc.ih;
      vec3<double> rh_calc;
      rh_calc = compute_h_position(rc, sites_cart);
      sites_cart[ih] = rh_calc;
    }
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
// INPUT: Gv, r1, r2
// Gr1, Gr2 are modified by reference
void G_crossproduct(
                   vec3<double> Gv,
                   vec3<double> r1,
                   vec3<double> r2,
                   vec3<double> & Gr1,
                   vec3<double> & Gr2)
{
  double x1 = r1[0], y1 = r1[1], z1 = r1[2];
  double x2 = r2[0], y2 = r2[1], z2 = r2[2];
  double Gv_x = Gv[0], Gv_y = Gv[1], Gv_z = Gv[2];
  Gr1[0] = -z2 * Gv_y + y2 * Gv_z;
  Gr1[1] =  z2 * Gv_x - x2 * Gv_z;
  Gr1[2] = -y2 * Gv_x + x2 * Gv_y;
  Gr2[0] =  z1 * Gv_y - y1 * Gv_z;
  Gr2[1] = -z1 * Gv_x + x1 * Gv_z;
  Gr2[2] =  y1 * Gv_x - x1 * Gv_y;
}

af::shared<scitbx::vec3<double> > modify_gradients_cpp(
                   af::shared<scitbx::vec3<double> > const& gradients,
                   af::shared<vec3<double> > const& sites_cart,
                   boost::python::list const& parameterization_)
  {
    af::shared<scitbx::vec3<double> >gradients_new(gradients.size());
    for(std::size_t i=0; i < gradients.size(); i++) {
      gradients_new[i] = gradients[i];
    }

    af::shared<riding_coefficients> parameterization;
    for(std::size_t i=0;i<boost::python::len(parameterization_);i++) {
      parameterization.push_back(
          boost::python::extract<riding_coefficients>(parameterization_[i])());
    }

    for(int i=0; i < parameterization.size(); i++) {
      riding_coefficients rc = parameterization[i];
      int ih = rc.ih, a0 = rc.a0, a1 = rc.a1, a2 = rc.a2, a3 = rc.a3;
      double a = rc.a, b = rc.b, h = rc.h, dh = rc.disth;
      vec3<double> r0 = sites_cart[a0];
      vec3<double> rh = sites_cart[rc.ih];
      vec3<double> r1 = sites_cart[rc.a1];
      vec3<double> GH = gradients_new[ih];
      vec3<double> G0 = gradients_new[a0];
      vec3<double> G1 = gradients_new[a1];
      vec3<double> G2 = gradients_new[a2];
      vec3<double> Gu10, Gu20, Gu30, G01, GuH0;
      vec3<double> Gr10, Gr20, Gr30;

      if (rc.htype == "flat_2neigbs" || rc.htype=="2neigbs" || rc.htype=="2tetra") {
        vec3<double> r2 = sites_cart[rc.a2];
        vec3<double> r10 = (r1-r0);
        vec3<double> r20 = (r2-r0);
        vec3<double> u10 = r10.normalize();
        vec3<double> u20 = r20.normalize();
        // step 1
        GuH0 = dh * GH;
        G01 = GH;
        // step 2
        if (rc.htype == "flat_2neigbs") {
          vec3<double> rh0 = (a * u10 + b * u20);
          double length = rh0.length();
          MMTBX_ASSERT(length > 0.);
          vec3<double>  multiplier1 = 1/length * GuH0;
          double multiplier2 = a * b * (pow(1/length,3)*(rh0*GuH0));
          Gu10 = a * multiplier1 - u20 * multiplier2;
          Gu20 = b * multiplier1 - u10 * multiplier2;
        } else if (rc.htype == "2tetra") {
          double alpha = h;
          vec3<double> v = u10.cross(u20);
          vec3<double> d = a * u10 + b * u20;
          // step 2
          vec3<double> Gd0 = cos(alpha) * GuH0;
          vec3<double> Gv0 = sin(alpha) * GuH0;
          // step 3
          vec3<double> Gd = G_unitvector(Gd0, d);
          vec3<double> Gv = G_unitvector(Gv0, v);
          // step 4
          vec3<double> Gu10_1 = a * Gd;
          vec3<double> Gu20_1 = b * Gd;
          // step 5
          vec3<double> Gu10_2, Gu20_2;
          G_crossproduct(Gv, u10, u20, Gu10_2, Gu20_2);
          Gu10 = Gu10_1 + Gu10_2;
          Gu20 = Gu20_1 + Gu20_2;
        } else if (rc.htype == "2neigbs") {
          vec3<double> v = u10.cross(u20);
          vec3<double> v0 = v.normalize();
          vec3<double> rh0 = a * u10 + b * u20 + h * v0;
          double length = rh0.length();
          MMTBX_ASSERT(length > 0.);
          vec3<double> multiplier1 = 1/length * GuH0;
          double multiplier2 = a * b* pow((1/length),3)*(rh0 * GuH0);
          vec3<double> Gu10_1 = a * multiplier1 - u20 * multiplier2;
          vec3<double> Gu20_1 = b * multiplier1 - u10 * multiplier2;
          vec3<double> Gv0 = h * multiplier1;
          // step 3
          vec3<double> Gv = G_unitvector(Gv0,v);
          // step 4
          vec3<double> Gu10_2, Gu20_2;
          G_crossproduct(Gv, u10, u20, Gu10_2, Gu20_2);
          Gu10 = Gu10_1 + Gu10_2;
          Gu20 = Gu20_1 + Gu20_2;
        }
        // step 5
        Gr10 = G_unitvector(Gu10, r10);
        Gr20 = G_unitvector(Gu20, r20);
        gradients_new[a0] = G0 + GH - Gr10 - Gr20;
        gradients_new[a1] = G1 + Gr10;
        gradients_new[a2] = G2 + Gr20;
      } else if (rc.htype == "3neigbs") {
        vec3<double> G3 = gradients_new[a3];
        vec3<double> r2 = sites_cart[a2];
        vec3<double> r3 = sites_cart[a3];
        vec3<double> r10 = (r1-r0), r20 = (r2-r0), r30 = (r3-r0);
        vec3<double> u10 = r10.normalize();
        vec3<double> u20 = r20.normalize();
        vec3<double> u30 = r30.normalize();
        vec3<double> rh0 = (a * u10 + b * u20 + h * u30);
        double length = rh0.length();
        MMTBX_ASSERT(length > 0.);
        // step 1
        GuH0 = dh * GH;
        G01 = GH;
        // step 2
        vec3<double> multiplier1 = 1/length * GuH0;
        double multiplier2 = pow((1/length),3) * (rh0*GuH0);
        Gu10 = a * multiplier1 - (a*b*u20+a*h*u30) * multiplier2;
        Gu20 = b * multiplier1 - (a*b*u10+b*h*u30) * multiplier2;
        Gu30 = h * multiplier1 - (a*h*u10+b*h*u20) * multiplier2;
        // step 3
        Gr10 = G_unitvector(Gu10, r10);
        Gr20 = G_unitvector(Gu20, r20);
        Gr30 = G_unitvector(Gu30, r30);
        gradients_new[a0] = G0 + GH - Gr10 - Gr20 -Gr30;
        gradients_new[a1] = G1 + Gr10;
        gradients_new[a2] = G2 + Gr20;
        gradients_new[a3] = G3 + Gr30;
      } else if (rc.htype == "alg1b" || rc.htype=="prop" || rc.htype=="alg1a") {
        vec3<double> rb1 = sites_cart[a2];
        double alpha = a, phi = b, n = rc.n;
        double k1 = -cos(alpha);
        double k2 = sin(alpha) * cos(phi + n * 2 * scitbx::constants::pi/3.0);
        double k3 = sin(alpha) * sin(phi + n * 2 * scitbx::constants::pi/3.0);
        vec3<double> r01 = r0 - r1, rb10 = rb1 - r1, v1 = r0 - r1;
        vec3<double> u1 = (r0 - r1).normalize();
        vec3<double> v2 = rb10 - (rb10 * u1) * u1;
        vec3<double> u2 = v2.normalize();
        vec3<double> u3 = u1.cross(u2);
        // step 1
        GuH0 = dh * GH;
        G01 = GH;
        // step 2
        vec3<double> Gu1_1 = k1 * GuH0;
        vec3<double> Gu2_1 = k2 * GuH0;
        vec3<double> Gu3   = k3 * GuH0;
        // step 3
        vec3<double> Gu1_2, Gu2_2;
        G_crossproduct(Gu3, u1, u2, Gu1_2, Gu2_2);
        vec3<double> Gu2 = Gu2_1 + Gu2_2;
        // step 4
        vec3<double> Gv2 = G_unitvector(Gu2, v2);
        // step 5
        vec3<double> Gu1 = Gu1_1 + Gu1_2 - ((rb10*u1) * Gv2) - ((u1*Gv2) * rb10);
        vec3<double> Grb1 = Gv2 - (u1 * Gv2) * u1;
        // step 6
        vec3<double> Gr01 = G_unitvector(Gu1, r01);
        gradients_new[a0] = G0 + GH + Gr01;
        gradients_new[a1] = G1 - Gr01 - Grb1;
        gradients_new[a2] = G2 + Grb1;
      }

    }
  return gradients_new;
} // modify_gradients_cpp


}} // namespace mmtbx::riding_h

#endif // MMTBX_RIDING_H
