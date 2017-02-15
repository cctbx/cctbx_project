#ifndef MMTBX_RIDING_H
#define MMTBX_RIDING_H

#include <iostream>
#include <string>
#include <cctbx/sgtbx/space_group.h>
#include <boost/python/list.hpp>

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

vec3<double> compute_H_position(
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
      double length = sqrt(a*a + b*b + 2*a*b*u10*u20);
      vec3<double> uh0 = (a * u10 + b * u20) / length;
      rh_calc = r0 + dh * uh0;
    } else if (rc.htype == "2neigbs") {
      vec3<double> r2 = sites_cart[rc.a2];
      vec3<double> u10 = (r1-r0).normalize();
      vec3<double> u20 = (r2-r0).normalize();
      vec3<double> v0 = (u10.cross(u20)).normalize();
      vec3<double> rh0 = (a * u10 + b * u20 + h * v0);
      double length = sqrt(rh0 * rh0);
      vec3<double> uh0 = rh0/length;
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
      double length = sqrt(rh0 * rh0);
      vec3<double> uh0 = rh0/length;
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
    std::cout << "length of parameterization: " << len(parameterization_) << "\n";
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
      //std::cout << "ih: " << ih << "  rx_ini: " << sites_cart[ih][0] << "\n";
      rh_calc = compute_H_position(rc, sites_cart);
      sites_cart_new[ih] = rh_calc;
    }
  return sites_cart_new;
}

//af::shared<scitbx::vec3<double> > modify_gradients_cpp(
//                   af::shared<scitbx::vec3<double> > const& gradients,
//                   boost::python::list const& parameterization_)
//  {
//    std::cout << "length of parameterization: " << len(parameterization_) << "\n";
//    af::shared<riding_coefficients> parameterization;
//    af::shared<scitbx::vec3<double> > gradients_new;
//    gradients_new = gradients;
//    for(std::size_t i=0;i<boost::python::len(parameterization_);i++) {
//      parameterization.push_back(
//          boost::python::extract<riding_coefficients>(parameterization_[i])());
//    }
//    for(int i=0; i < parameterization.size(); i++) {
//      riding_coefficients rc = parameterization[i];
//      int ih = rc.ih;
//      std::cout << "ih: " << ih << "  Gx: " << gradients[ih][0] << "\n";
//    }
//  return gradients_new;
//}


}} // namespace mmtbx::riding_h

#endif // MMTBX_RIDING_H
