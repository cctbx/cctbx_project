#ifndef MMTBX_TLS_H
#define MMTBX_TLS_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <vector>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/versa_algebra.h>
#include <cctbx/import_scitbx_af.h>
#include <cmath>
#include <cctbx/adptbx.h>
#include <cctbx/xray/scattering_type_registry.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <mmtbx/error.h>
#include <cctbx/xray/targets.h>
#include <scitbx/matrix/outer_product.h>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>

using namespace std;
namespace mmtbx { namespace tls {
namespace af=scitbx::af;
using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;

class d_uaniso_d_tls {
public:
  d_uaniso_d_tls(vec3<double> const& origin,
                 vec3<double> const& site,
                 bool scale_l_and_s,
                 bool use_trace_s_zero_constraint)
  {
    // Order is : sym_mat3<double>(a11,a22,a33,a12,a13,a23);
    //                            ( a0, a1, a2, a3, a4, a5);
    // For general 3x3 matrix: | 0 1 2 |
    //                         | 3 4 5 | = S
    //                         | 6 7 8 |
    if(scale_l_and_s) {
       s_sc = scitbx::deg_as_rad(1.0);
       l_sc = s_sc * s_sc;
    }
    double const nul=0.0, one=1.0;
    vec3<double> r = site - origin;
    x = r[0];
    y = r[1];
    z = r[2];
    xx = x*x;
    yy = y*y;
    zz = z*z;
    xy = x*y;
    yz = y*z;
    xz = x*z;
    xy2= -2.*xy;
    xz2= -2.*xz;
    yz2= -2.*yz;
    z2 = 2.*z;
    y2 = 2.*y;
    x2 = 2.*x;
    gT0 = sym_mat3<double>(one, nul, nul, nul, nul, nul);
    gT1 = sym_mat3<double>(nul, one, nul, nul, nul, nul);
    gT2 = sym_mat3<double>(nul, nul, one, nul, nul, nul);
    gT3 = sym_mat3<double>(nul, nul, nul, one, nul, nul);
    gT4 = sym_mat3<double>(nul, nul, nul, nul, one, nul);
    gT5 = sym_mat3<double>(nul, nul, nul, nul, nul, one);
    gL0 = sym_mat3<double>(nul,  zz,  yy, nul, nul, -yz);
    gL1 = sym_mat3<double>( zz, nul,  xx, nul, -xz, nul);
    gL2 = sym_mat3<double>( yy,  xx, nul, -xy, nul, nul);
    gL3 = sym_mat3<double>(nul, nul, xy2, -zz,  yz,  xz);
    gL4 = sym_mat3<double>(nul, xz2, nul,  yz, -yy,  xy);
    gL5 = sym_mat3<double>(yz2, nul, nul,  xz,  xy, -xx);
    if(use_trace_s_zero_constraint) {
       gS0 = sym_mat3<double>(nul, nul, nul,  -z,  y2, -x);
       gS4 = sym_mat3<double>(nul, nul, nul,   z,   y, -x2);
       gS8 = sym_mat3<double>(nul, nul, nul, nul, nul, nul);
    }
    else {
       gS0 = sym_mat3<double>(nul, nul, nul,  -z,   y, nul);
       gS4 = sym_mat3<double>(nul, nul, nul,   z, nul,  -x);
       gS8 = sym_mat3<double>(nul, nul, nul, nul,  -y,   x);
    }
    gS1 = sym_mat3<double>(nul, -z2, nul, nul, nul,   y);
    gS2 = sym_mat3<double>(nul, nul,  y2, nul, nul,  -z);
    gS3 = sym_mat3<double>( z2, nul, nul, nul,  -x, nul);
    gS5 = sym_mat3<double>(nul, nul, -x2, nul,   z, nul);
    gS6 = sym_mat3<double>(-y2, nul, nul,   x, nul, nul);
    gS7 = sym_mat3<double>(nul,  x2, nul,  -y, nul, nul);
    gT_ = sym_mat3<sym_mat3<double> >(gT0,gT1,gT2,gT3,gT4,gT5);
    if(scale_l_and_s) {
       gL_ = sym_mat3<sym_mat3<double> >
                       (gL0*l_sc,gL1*l_sc,gL2*l_sc,gL3*l_sc,gL4*l_sc,gL5*l_sc);
       gS_ = mat3<sym_mat3<double> >
            (gS0*s_sc,gS1*s_sc,gS2*s_sc,gS3*s_sc,gS4*s_sc,gS5*s_sc,gS6*s_sc,
                                                            gS7*s_sc,gS8*s_sc);
    }
    else {
       gL_ = sym_mat3<sym_mat3<double> >(gL0,gL1,gL2,gL3,gL4,gL5);
       gS_ = mat3<sym_mat3<double> >(gS0,gS1,gS2,gS3,gS4,gS5,gS6,gS7,gS8);
    }
  }
  sym_mat3<sym_mat3<double> > d_u_d_T() { return gT_; }
  sym_mat3<sym_mat3<double> > d_u_d_L() { return gL_; }
  mat3<sym_mat3<double> >     d_u_d_S() { return gS_; }
private:
  double x,y,z,xx,yy,zz,xy,yz,xz,xy2,xz2,yz2,z2,y2,x2,s_sc,l_sc;
  sym_mat3<double> gT0,gT1,gT2,gT3,gT4,gT5;
  sym_mat3<double> gL0,gL1,gL2,gL3,gL4,gL5;
  sym_mat3<double> gS0,gS1,gS2,gS3,gS4,gS5,gS6,gS7,gS8;
  sym_mat3<sym_mat3<double> > gT_,gL_;
  mat3<sym_mat3<double> > gS_;
};

sym_mat3<double> u_cart_from_xyz(af::shared<vec3<double> > const& sites_cart)
{
  int size = sites_cart.size();
  double cmx=0;
  double cmy=0;
  double cmz=0;
  for(std::size_t i=0; i < size; i++) {
    cmx += sites_cart[i][0];
    cmy += sites_cart[i][1];
    cmz += sites_cart[i][2];
  }
  vec3<double> cm(cmx/size, cmy/size, cmz/size);
  af::shared<double> x(size);
  af::shared<double> y(size);
  af::shared<double> z(size);
  af::shared<vec3<double> > sites_cart_cm = sites_cart-cm;
  for(std::size_t i=0; i < sites_cart_cm.size(); i++) {
    x[i]=sites_cart_cm[i][0];
    y[i]=sites_cart_cm[i][1];
    z[i]=sites_cart_cm[i][2];
  }
  double u11 = 0.;
  double u22 = 0.;
  double u33 = 0.;
  double u12 = 0.;
  double u13 = 0.;
  double u23 = 0.;
  for(std::size_t i=0; i < sites_cart_cm.size(); i++) {
    u11 += x[i]*x[i];
    u22 += y[i]*y[i];
    u33 += z[i]*z[i];
    u12 += x[i]*y[i];
    u13 += x[i]*z[i];
    u23 += y[i]*z[i];
  }
  return sym_mat3<double>(u11/size,u22/size,u33/size,u12/size,u13/size,u23/size);;
}

boost::python::list all_vs_all(
  boost::python::list const& xyz_all_
  )

{
  boost::python::list xyz_atoms_all;
  af::shared< af::shared<vec3<double> > > xyz_all;

  for(std::size_t i=0;i<boost::python::len(xyz_all_);i++) {
    xyz_all.push_back(
        boost::python::extract<af::shared<vec3<double> > >(xyz_all_[i])());
  }


  //af::shared< af::shared<vec3<double> > > xyz_atoms_all(xyz_all[0].size());
  for(std::size_t i=0; i < xyz_all[0].size(); i++) {
    af::shared<vec3<double> > xyz_atoms;
    for(std::size_t j=0; j < xyz_all.size(); j++) {
      xyz_atoms.push_back(xyz_all[j][i]);
    }
    //xyz_atoms_all.push_back(xyz_atoms);
    xyz_atoms_all.append(xyz_atoms);
  }

  return xyz_atoms_all;
}


af::shared<vec3<double> > apply_tls_shifts(
  af::shared<vec3<double> > const& sites_cart,
  mat3<double> const& R_ML_transposed,
  mat3<double> const& R_ML,
  vec3<double> const& d0,
  vec3<double> const& d_r_M_V,
  vec3<double> const& s,
  double const& wy_lx,
  double const& wz_lx,
  double const& wz_ly,
  double const& wx_ly,
  double const& wx_lz,
  double const& wy_lz,
  vec3<double> const& origin
  )
{
  af::shared<vec3<double> > sites_cart_new(sites_cart.size());
  for(std::size_t i=0; i < sites_cart.size(); i++) {
    vec3<double> r_L = R_ML_transposed * sites_cart[i];
    //
    double sxb = s[0];
    double syb = s[1];
    double szb = s[2];
    double dx0 = d0[0];
    double dy0 = d0[1];
    double dz0 = d0[2];
    double x = r_L[0];
    double y = r_L[1];
    double z = r_L[2];
    vec3<double> d_lx_r_L = vec3<double>(
      sxb*dx0,
      (y-wy_lx)*(std::cos(dx0)-1) - (z-wz_lx)*std::sin(dx0),
      (y-wy_lx)*std::sin(dx0)     + (z-wz_lx)*(std::cos(dx0)-1)
      );
    vec3<double> d_ly_r_L = vec3<double> (
      (z-wz_ly)*std::sin(dy0)     + (x-wx_ly)*(std::cos(dy0)-1),
      syb*dy0,
      (z-wz_ly)*(std::cos(dy0)-1) - (x-wx_ly)*std::sin(dy0)
      );
    vec3<double> d_lz_r_L = vec3<double> (
      (x-wx_lz)*(std::cos(dz0)-1) - (y-wy_lz)*std::sin(dz0),
      (x-wx_lz)*std::sin(dz0)     + (y-wy_lz)*(std::cos(dz0)-1),
      szb*dz0
      );
    vec3<double> d_r_L = d_lx_r_L + d_ly_r_L + d_lz_r_L;
    vec3<double> d_r_M_L = R_ML * d_r_L;
    //
    vec3<double> d_r_M = d_r_M_L + d_r_M_V;
    sites_cart_new[i] = sites_cart[i] + d_r_M + origin;
  }

  return sites_cart_new;
}

double uiso_from_tls(double const& T,
                     sym_mat3<double> const& L_deg,
                     vec3<double> const& S_deg,
                     vec3<double> const& origin,
                     vec3<double> const& site_cart)
  {
   double deg2rad = scitbx::deg_as_rad(1.0);
   double deg2radsq = deg2rad * deg2rad;
   sym_mat3<double> L = L_deg * deg2radsq;
   vec3<double> S = S_deg * deg2rad;
   vec3<double> r = site_cart - origin;
   double x = r[0];
   double y = r[1];
   double z = r[2];
   double u_iso = T + 1./3*(L[0]*(y*y+z*z) + L[1]*(x*x+z*z) + L[2]*(x*x+y*y) -
                            2*L[3]*x*y - 2*L[4]*x*z - 2*L[5]*y*z +
                            2*S[0]*z   + 2*S[1]*y   + 2*S[2]*x);
  return u_iso;
}

class d_uiso_d_tls {
public:
  d_uiso_d_tls(vec3<double> const& origin,
               vec3<double> const& site,
               double scale_l_and_s)
  {
   // Order is : sym_mat3<double>(a11,a22,a33,a12,a13,a23);
   //                            ( a0, a1, a2, a3, a4, a5);
   vec3<double> r = site - origin;
   double x = r[0];
   double y = r[1];
   double z = r[2];
   double one_over_three = 1./3;
   double two_over_three = 2./3;
   d_uiso_d_T = 1;
   double d_uiso_d_L11 = one_over_three * (y*y+z*z);
   double d_uiso_d_L22 = one_over_three * (x*x+z*z);
   double d_uiso_d_L33 = one_over_three * (x*x+y*y);
   double d_uiso_d_L12 = -two_over_three * x*y;
   double d_uiso_d_L13 = -two_over_three * x*z;
   double d_uiso_d_L23 = -two_over_three * y*z;
   double d_uiso_d_S1  =  two_over_three * z;
   double d_uiso_d_S2  =  two_over_three * y;
   double d_uiso_d_S3  =  two_over_three * x;
   if(scale_l_and_s != -1) {
     double l_sc = scale_l_and_s * scale_l_and_s;
     d_uiso_d_L = sym_mat3<double>
       (d_uiso_d_L11*l_sc,
        d_uiso_d_L22*l_sc,
        d_uiso_d_L33*l_sc,
        d_uiso_d_L12*l_sc,
        d_uiso_d_L13*l_sc,
        d_uiso_d_L23*l_sc);
     d_uiso_d_S = vec3<double> (
       d_uiso_d_S1*scale_l_and_s,
       d_uiso_d_S2*scale_l_and_s,
       d_uiso_d_S3*scale_l_and_s);
   }
   else {
    d_uiso_d_L = sym_mat3<double>
       (d_uiso_d_L11,
        d_uiso_d_L22,
        d_uiso_d_L33,
        d_uiso_d_L12,
        d_uiso_d_L13,
        d_uiso_d_L23);
    d_uiso_d_S = vec3<double> (d_uiso_d_S1, d_uiso_d_S2, d_uiso_d_S3);
   }
  }
  double d_u_d_T() { return d_uiso_d_T; }
  sym_mat3<double> d_u_d_L() { return d_uiso_d_L; }
  vec3<double> d_u_d_S() { return d_uiso_d_S; }
private:
  double d_uiso_d_T;
  sym_mat3<double> d_uiso_d_L;
  vec3<double> d_uiso_d_S;
};

class tls_from_uiso_target_and_grads {
public:
  tls_from_uiso_target_and_grads(
                                  double const& T,
                                  sym_mat3<double> const& L_deg,
                                  vec3<double> const& S_deg,
                                  vec3<double> const& origin,
                                  af::shared<vec3<double> > const& sites,
                                  af::shared<double> const& uisos)
  {
    tg = 0.0;
    gT = 0.0;
    gL = sym_mat3<double> (0,0,0,0,0,0);
    gS = vec3<double> (0,0,0);
    double scale_l_and_s = scitbx::deg_as_rad(1.0);
    for(std::size_t i=0; i < sites.size(); i++) {
      vec3<double> const& site = sites[i];
      double diff = uiso_from_tls(T,L_deg,S_deg,origin,site) - uisos[i];
      double two_diff = 2*diff;
      tg += diff*diff;
      d_uiso_d_tls d_uiso_d_tls_manager(origin,site,scale_l_and_s);
      gT = gT + two_diff*d_uiso_d_tls_manager.d_u_d_T();
      gL = gL + two_diff*d_uiso_d_tls_manager.d_u_d_L();
      gS = gS + two_diff*d_uiso_d_tls_manager.d_u_d_S();
    }
  }
  double target() const { return tg; }
  double grad_T() { return gT; }
  sym_mat3<double> grad_L() { return gL; }
  vec3<double> grad_S() { return gS; }
private:
  double tg, gT;
  sym_mat3<double> gL;
  vec3<double> gS;
};

class uaniso_from_tls {
public:
  uaniso_from_tls(sym_mat3<double> const& T,
                  sym_mat3<double> const& L_deg,
                  mat3<double> const& S_deg,
                  vec3<double> const& origin,
                  vec3<double> const& site_cart,
                  bool zeroize_trace)
  {
   double deg2rad   = scitbx::deg_as_rad(1.0);
   double deg2radsq = deg2rad * deg2rad;
   L = L_deg * deg2radsq;
   S = S_deg * deg2rad;
   if(zeroize_trace) S[8] = -(S[0]+S[4]); // condition on trace(S) = 0.0
   r_ = site_cart - origin;
   x = r_[0];
   y = r_[1];
   z = r_[2];
   ALA_ = L.antisymmetric_tensor_transform(z,-y,x);
   ASSA_ = sym_mat3<double>(2.*S[3]*z - 2.*S[6]*y,
                            2.*S[7]*x - 2.*S[1]*z,
                            2.*S[2]*y - 2.*S[5]*x,
                            S[6]*x - S[7]*y + (S[4] - S[0])*z,
                            S[5]*z - S[3]*x + (S[0] - S[8])*y,
                            S[1]*y - S[2]*z + (S[8] - S[4])*x);
   uaniso = sym_mat3<double>(T + ALA_ + ASSA_);
  }
  double x,y,z;
  sym_mat3<double> u() { return uaniso; }
  sym_mat3<double> Lrad() { return L; }
  mat3<double> Srad() { return S; }
  vec3<double> r() { return r_; }
  sym_mat3<double> ALA() { return ALA_; }
  sym_mat3<double> ASSA() { return ASSA_; }
private:
  sym_mat3<double> uaniso;
  sym_mat3<double> L;
  mat3<double> S;
  vec3<double> r_;
  sym_mat3<double> ALA_;
  sym_mat3<double> ASSA_;
};

//
class d_target_d_tls {
public:
  d_target_d_tls(af::shared<vec3<double> > const& sites,
                 vec3<double> const& origin,
                 af::shared<sym_mat3<double> > const& d_target_d_uaniso,
                 bool scale_l_and_s,
                 bool use_trace_s_zero_constraint)
  {
    gT.resize(6,0.0);
    gL.resize(6,0.0);
    gS.resize(9,0.0);
    for (std::size_t i=0; i < sites.size(); i++) {
         vec3<double> const& site = sites[i];
         d_uaniso_d_tls d_uaniso_d_tls_manager(origin,site,scale_l_and_s,
                                                  use_trace_s_zero_constraint);
         sym_mat3<sym_mat3<double> > d_u_d_T =d_uaniso_d_tls_manager.d_u_d_T();
         sym_mat3<sym_mat3<double> > d_u_d_L =d_uaniso_d_tls_manager.d_u_d_L();
         mat3<sym_mat3<double> >     d_u_d_S =d_uaniso_d_tls_manager.d_u_d_S();
         for (std::size_t k=0; k < 6; k++) {
              for (std::size_t m=0; m < 6; m++) {
                   if(d_u_d_T[k][m] != 0.0) {
                      gT[k] += d_target_d_uaniso[i][m] * d_u_d_T[k][m];
                   }
                   if(d_u_d_L[k][m] != 0.0) {
                      gL[k] += d_target_d_uaniso[i][m] * d_u_d_L[k][m];
                   }
                   if(d_u_d_S[k][m] != 0.0) {
                      gS[k] += d_target_d_uaniso[i][m] * d_u_d_S[k][m];
                   }
             }
         }
         for (std::size_t k=6; k < 9; k++) {
              for (std::size_t m=0; m < 6; m++) {
                   if(d_u_d_S[k][m] != 0.0) {
                      gS[k] += d_target_d_uaniso[i][m] * d_u_d_S[k][m];
                   }
             }
         }
    }
  }
  af::shared<double> grad_T() { return gT; }
  af::shared<double> grad_L() { return gL; }
  af::shared<double> grad_S() { return gS; }
private:
  af::shared<double> gT, gL, gS;
};

class tls_from_uaniso_target_and_grads {
public:
  tls_from_uaniso_target_and_grads(
                                  sym_mat3<double> const& T,
                                  sym_mat3<double> const& L_deg,
                                  mat3<double> const& S_deg,
                                  vec3<double> const& origin,
                                  af::shared<vec3<double> > const& sites,
                                  af::shared<sym_mat3<double> > const& uanisos)
  {
    tg = 0.0;
    for (std::size_t i=0; i < sites.size(); i++) {
         vec3<double> const& site = sites[i];
         uaniso_from_tls manager(T,L_deg,S_deg,origin,site,true);
         sym_mat3<double> const& utls = manager.u();
         sym_mat3<double> diff = utls - uanisos[i];
         for (std::size_t k=0; k < diff.size(); k++) {
              tg += diff[k] * diff[k];
         }
         diffs.push_back(diff*2.);
    }
    d_target_d_tls  d_target_d_tls_manager(sites, origin, diffs, true, false);
    gT = d_target_d_tls_manager.grad_T();
    gL = d_target_d_tls_manager.grad_L();
    gS = d_target_d_tls_manager.grad_S();
  }
  double target() const { return tg; }
  af::shared<double> grad_T() { return gT; }
  af::shared<double> grad_L() { return gL; }
  af::shared<double> grad_S() { return gS; }
private:
  double tg;
  af::shared<double> gT, gL, gS;
  af::shared<sym_mat3<double> > diffs;
};

template <typename FloatType=double>
class tlso
{
  public:
    scitbx::sym_mat3<FloatType> t;
    scitbx::sym_mat3<FloatType> l;
    scitbx::mat3<FloatType> s;
    scitbx::vec3<FloatType> origin;
    tlso() {
     t.fill(0.0);
     l.fill(0.0);
     s.fill(0.0);
     origin.fill(0.0);
    }
    tlso(scitbx::sym_mat3<FloatType> const& t_,
         scitbx::sym_mat3<FloatType> const& l_,
         scitbx::mat3<FloatType> const& s_,
         scitbx::vec3<FloatType> const&  origin_)
    :
      t(t_),l(l_),s(s_),origin(origin_)
    {}
};

class tls_parts_one_group {
public:
  tls_parts_one_group(tlso<double> tls_params,
                      af::shared<vec3<double> > const& sites_cart)
  {
   for (std::size_t i=0; i < sites_cart.size(); i++) {
        vec3<double> const& site = sites_cart[i];
        uaniso_from_tls manager(tls_params.t,
                                tls_params.l,
                                tls_params.s,
                                tls_params.origin,
                                site,
                                true);
        ala_.push_back(manager.ALA());
        assa_.push_back(manager.ASSA());
        u_cart_.push_back(manager.u());
        r_.push_back(manager.r());
        t_.push_back(tls_params.t);

   }
 }
  af::shared<sym_mat3<double> > ala() { return ala_; }
  af::shared<sym_mat3<double> > t() { return t_; }
  af::shared<sym_mat3<double> > assa() { return assa_; }
  af::shared<sym_mat3<double> > u_cart() { return u_cart_; }
  af::shared<vec3<double> > r() { return r_; }
protected:
  af::shared<sym_mat3<double> > ala_;
  af::shared<sym_mat3<double> > t_;
  af::shared<sym_mat3<double> > assa_;
  af::shared<sym_mat3<double> > u_cart_;
  af::shared<vec3<double> > r_;
};

class tls_parts_one_group_as_b_iso {
public:
  tls_parts_one_group_as_b_iso(tlso<double> tls_params,
                               af::shared<vec3<double> > const& sites_cart)
  {
   for (std::size_t i=0; i < sites_cart.size(); i++) {
        vec3<double> const& site = sites_cart[i];
        uaniso_from_tls manager(tls_params.t,
                                tls_params.l,
                                tls_params.s,
                                tls_params.origin,
                                site,
                                true);
        ala_.push_back(   cctbx::adptbx::u_as_b(cctbx::adptbx::u_cart_as_u_iso(manager.ALA())));
        assa_.push_back(  cctbx::adptbx::u_as_b(cctbx::adptbx::u_cart_as_u_iso(manager.ASSA())));
        b_iso_.push_back(cctbx::adptbx::u_as_b(cctbx::adptbx::u_cart_as_u_iso(manager.u())));
        t_.push_back(     cctbx::adptbx::u_as_b(cctbx::adptbx::u_cart_as_u_iso(tls_params.t)));
   }
 }
  af::shared<double> ala() { return ala_; }
  af::shared<double> t() { return t_; }
  af::shared<double> assa() { return assa_; }
  af::shared<double> b_iso() { return b_iso_; }
protected:
  af::shared<double> ala_;
  af::shared<double> t_;
  af::shared<double> assa_;
  af::shared<double> b_iso_;
};

class common {
  public:
    sym_mat3<double> t_;
    common() {}
    common(sym_mat3<double> const& m_,
           sym_mat3<double> const& n_,
           double small = 1.e-9)
    :
      t_(-1,-1,-1,-1,-1,-1),
      t1(-1,-1,-1,-1,-1,-1),
      t2(-1,-1,-1,-1,-1,-1),
      t3(-1,-1,-1,-1,-1,-1),
      g(-1,-1,-1,-1,-1,-1,-1,-1,-1),
      p(-1,-1,-1,-1,-1,-1,-1,-1,-1),
      m(m_),m_cp(m_),
      n(n_),n_cp(n_),
      small(small),small6(1.e-6),small5(1.e-5),small4(1.e-4),small3(1.e-1),
      a(0),bn(0),
      branch_0(false),branch_1(false),branch_1_1(false),branch_1_2(false),
      branch_1_2_1(false),branch_1_2_2(false),branch_1_2_3(false),
      branch_1_2_3_1(false),branch_1_2_3_2(false)
    {
      func_1();
      MMTBX_ASSERT((branch_0 || branch_1) && (branch_0 != branch_1));
      if(branch_1) {
         bool not_sure_2 = func_2();
         MMTBX_ASSERT((branch_1_1 || branch_1_2)&&(branch_1_1 != branch_1_2));
         if(not_sure_2 == false) {
            if(branch_1_1) {
               func_branch_1_1();
            }
            if(branch_1_2) {
               func_branch_1_2();
            }
         }
         else {
            if(branch_1_1) {
               try {
                 func_branch_1_1();
               }
               catch(const error& e) {
                 re_init_all();
                 func_branch_1_2();
               }
            }
            else {
               try {
                 func_branch_1_2();
               }
               catch(const error& e) {
                 re_init_all();
                 func_branch_1_1();
               }
            }
         }
      }

      bool test1 = cctbx::adptbx::is_positive_definite(t_,      small3);
      bool test2 = cctbx::adptbx::is_positive_definite(m_cp-t_, small3);
      bool test3 = cctbx::adptbx::is_positive_definite(n_cp-t_, small3);
      if(!test1 || !test2 || !test3) show_all();
      MMTBX_ASSERT( test1 );
      MMTBX_ASSERT( test2 );
      MMTBX_ASSERT( test3 );
    }

    void func_branch_1_1()
    {
      func_3();
      func_4();
      func_12();
    }

    void func_branch_1_2()
    {
      func_5();
      MMTBX_ASSERT(branch_1_2_1 || branch_1_2_2 || branch_1_2_3);
      if(branch_1_2_1) {
         MMTBX_ASSERT(!branch_1_2_2 && !branch_1_2_3);
         func_6();
         func_11();
         func_12();
      }
      if(branch_1_2_2) {
         MMTBX_ASSERT(!branch_1_2_1 && !branch_1_2_3);
         func_7();
         func_8();
         func_12();
      }
      if(branch_1_2_3) {
         MMTBX_ASSERT(!branch_1_2_1 && !branch_1_2_2);
         func_9();
         MMTBX_ASSERT((branch_1_2_3_1 || branch_1_2_3_2) &&
                      (branch_1_2_3_1 != branch_1_2_3_2));
         if(branch_1_2_3_1) {
            func_8();
            func_11();
            func_12();
         }
         if(branch_1_2_3_2) {
            func_10();
            func_4();
            func_12();
         }
      }
    }

    void re_init_all()
    {
      t_ = sym_mat3<double> (-1,-1,-1,-1,-1,-1);
      t1 = sym_mat3<double> (-1,-1,-1,-1,-1,-1);
      t2 = sym_mat3<double> (-1,-1,-1,-1,-1,-1);
      g = mat3<double> (-1,-1,-1,-1,-1,-1,-1,-1,-1),
      p = mat3<double> (-1,-1,-1,-1,-1,-1,-1,-1,-1),
      m = m_cp-t3;
      n = n_cp-t3;
      a = 0.0;
      bn = 0.0;
      branch_0=false, branch_1=false, branch_1_1=false, branch_1_2=false;
      branch_1_2_1=false, branch_1_2_2=false, branch_1_2_3=false;
      branch_1_2_3_1=false, branch_1_2_3_2=false;
    }

    void func_1()
     {
       cctbx::adptbx::eigensystem<double> m_eigensystem(m);
       cctbx::adptbx::eigensystem<double> n_eigensystem(n);
       vec3<double> myu = m_eigensystem.values();
       vec3<double> nyu = n_eigensystem.values();
       MMTBX_ASSERT( cctbx::adptbx::is_positive_definite(myu, small) );
       MMTBX_ASSERT( cctbx::adptbx::is_positive_definite(nyu, small) );
       MMTBX_ASSERT(myu[0]>=myu[1] && myu[1]>=myu[2] && myu[2]>0.);
       MMTBX_ASSERT(nyu[0]>=nyu[1] && nyu[1]>=nyu[2] && nyu[2]>0.);
       if(myu[2]<nyu[2]||(std::abs(myu[2]-nyu[2])<small && myu[1]<nyu[1])||
          (std::abs(myu[2]-nyu[2])<small&&std::abs(myu[1]-nyu[1])<small &&
          myu[0]<nyu[0])) {
          sym_mat3<double> tmp = m;
          m = n;
          n = tmp;
       }
       cctbx::adptbx::eigensystem<double> m_eigensystem_(m);
       cctbx::adptbx::eigensystem<double> n_eigensystem_(n);
       myu = m_eigensystem_.values();
       nyu = n_eigensystem_.values();
       if(!(myu[2] < nyu[0])) {
         t_ = n;
         branch_0 = true;
       }
       else branch_1 = true;
     }

     bool func_2()
      {
        cctbx::adptbx::eigensystem<double> m_eigensystem_(m);
        cctbx::adptbx::eigensystem<double> n_eigensystem_(n);
        vec3<double> myu = m_eigensystem_.values();
        vec3<double> nyu = n_eigensystem_.values();
        MMTBX_ASSERT(branch_1 != false);
        t3 = sym_mat3<double>(nyu[2],nyu[2],nyu[2],0,0,0);
        m = m - t3;
        n = n - t3;
        cctbx::adptbx::eigensystem<double> m_eigensystem(m);
        cctbx::adptbx::eigensystem<double> n_eigensystem(n);
        myu = m_eigensystem.values();
        nyu = n_eigensystem.values();
        vec3<double> g1 = n_eigensystem.vectors(0);
        vec3<double> g2 = n_eigensystem.vectors(1);
        vec3<double> g3 = n_eigensystem.vectors(2);
        double tmp = g3 * (m * g3);
        if(tmp > 0.0 && std::abs(tmp) < small) tmp = 0.0;
        double tmp2 = std::abs(nyu[1]);
        if(std::abs(tmp2) < small) tmp2 = 0.0;
        if(tmp2 > 0 && std::abs(myu[2]) < small && tmp > 0) {
           branch_1_1 = true;
        }
        else {
           branch_1_2 = true;
        }
        bool not_sure = false;
        if(std::abs(tmp) < small4 && std::abs(myu[2]) < small4 &&
           std::abs(tmp) > small  && std::abs(myu[2]) > small ) not_sure=true;
        return not_sure;
      }

     void func_3()
     {
       g = mat3<double>(1,0,0, 0,1,0, 0,0,1);
       cctbx::adptbx::eigensystem<double> m_eigensystem(m);
       cctbx::adptbx::eigensystem<double> n_eigensystem(n);
       vec3<double> myu = m_eigensystem.values();
       vec3<double> nyu = n_eigensystem.values();
       MMTBX_ASSERT(myu[0]>=myu[1]&&myu[1]>=myu[2]&&std::abs(myu[2])<small4);
       MMTBX_ASSERT(nyu[0]>=nyu[1]&&nyu[1]> nyu[2]&&std::abs(nyu[2])<small4);
       vec3<double> e3 = m_eigensystem.vectors(2);
       vec3<double> g3 = n_eigensystem.vectors(2);
       MMTBX_ASSERT(std::abs(g3[0]-e3[0])>small || std::abs(g3[1]-e3[1])>small
                    || std::abs(g3[2]-e3[2])>small);
       a = 0.0;
       t2 = sym_mat3<double>(0,0,0,0,0,0);
     }

     void func_4()
     {
       cctbx::adptbx::eigensystem<double> m_eigensystem(m);
       vec3<double> e3 = m_eigensystem.vectors(2);
       cctbx::adptbx::eigensystem<double> n_eigensystem(n);
       vec3<double> g3 = n_eigensystem.vectors(2);
       vec3<double> p3 = g3;
       vec3<double> u  = e3-(e3*p3)*p3;
       vec3<double> p2 = u / u.length();
       vec3<double> p1 = p2.cross(p3);
       double p1_n = p1.length();
       double p2_n = p2.length();
       double p3_n = p3.length();
       p = mat3<double>(p1[0]/p1_n,p2[0]/p2_n,p3[0]/p3_n,
                        p1[1]/p1_n,p2[1]/p2_n,p3[1]/p3_n,
                        p1[2]/p1_n,p2[2]/p2_n,p3[2]/p3_n);
       m = sym_mat3<double>(p.transpose() * m * p, small);
       n = sym_mat3<double>(p.transpose() * n * p, small);
       bn = beta(n, small);
     }

     void func_5()
     {
       cctbx::adptbx::eigensystem<double> n_eigensystem(n);
       vec3<double> g1 = n_eigensystem.vectors(0);
       vec3<double> g2 = n_eigensystem.vectors(1);
       vec3<double> g3 = n_eigensystem.vectors(2);
       g = mat3<double>(g1[0],g2[0],g3[0],
                        g1[1],g2[1],g3[1],
                        g1[2],g2[2],g3[2]);
       vec3<double> nyu = n_eigensystem.values();
       n = sym_mat3<double>(nyu[0],nyu[1],0,0,0,0);
       m = sym_mat3<double>(g.transpose() * m * g, small);
       cctbx::adptbx::eigensystem<double> m_eigensystem(m);
       cctbx::adptbx::eigensystem<double> n_eigensystem_(n);
       vec3<double> myu = m_eigensystem.values();
       nyu = n_eigensystem_.values();
       if(std::abs(nyu[1])<small) branch_1_2_1 = true;
       else {
          if(std::abs(myu[2])<small) branch_1_2_2 = true;
          else branch_1_2_3 = true;
       }
     }

     void func_6()
     {
       MMTBX_ASSERT(branch_1_2_1 != false);
       cctbx::adptbx::eigensystem<double> m_eigensystem(m);
       cctbx::adptbx::eigensystem<double> n_eigensystem(n);
       vec3<double> myu = m_eigensystem.values();
       vec3<double> nyu = n_eigensystem.values();
       MMTBX_ASSERT(myu[0]>=myu[1]&&myu[1]>=myu[2]&&(std::abs(myu[2])<small ||
                    myu[2]>0.0));
       MMTBX_ASSERT(nyu[0]>=nyu[1]&&std::abs(nyu[1])<small&&
                    std::abs(nyu[2])<small);
       a = 0.0;
       t2 = sym_mat3<double>(0,0,0,0,0,0);
     }

     void func_7()
     {
       MMTBX_ASSERT(branch_1_2_2 != false);
       cctbx::adptbx::eigensystem<double> m_eigensystem(m);
       cctbx::adptbx::eigensystem<double> n_eigensystem(n);
       vec3<double> myu = m_eigensystem.values();
       vec3<double> nyu = n_eigensystem.values();
       MMTBX_ASSERT(myu[0]>=myu[1]&&myu[1]>=myu[2]&&std::abs(myu[2])<small);
       MMTBX_ASSERT(nyu[0]>=nyu[1]&&nyu[1]> nyu[2]&&std::abs(nyu[2])<small);
       vec3<double> e3 = m_eigensystem.vectors(2);
       vec3<double> g3 = n_eigensystem.vectors(2);
       MMTBX_ASSERT(std::abs(g3[0]-e3[0])<small3&&std::abs(g3[1]-e3[1])<small3
                    && std::abs(g3[2]-e3[2])<small3);
       a = nyu[1];
     }

     void func_8()
     {
       t2 = sym_mat3<double>(a,a,0,0,0,0);
       cctbx::adptbx::eigensystem<double> n_eigensystem(n);
       vec3<double> nyu = n_eigensystem.values();
       n = sym_mat3<double>(nyu[0]-nyu[1],0,0,0,0,0);
       m = m - t2;
     }

     void func_9()
     {
       cctbx::adptbx::eigensystem<double> m_eigensystem(m);
       cctbx::adptbx::eigensystem<double> n_eigensystem(n);
       vec3<double> myu = m_eigensystem.values();
       vec3<double> nyu = n_eigensystem.values();
       MMTBX_ASSERT(myu[0]>=myu[1]&&myu[1]>=myu[2]&&myu[2]>0.0);
       MMTBX_ASSERT(nyu[0]>=nyu[1]&&nyu[1]> nyu[2]&&std::abs(nyu[2])<small);
       double am = alpha(m);
       if(nyu[1] <= am) {
          a = nyu[1];
          branch_1_2_3_1 = true;
       }
       else {
          a = am;
          branch_1_2_3_2 = true;
       }
     }

     void func_10()
     {
       t2 = sym_mat3<double>(a,a,0,0,0,0);
       m = m - t2;
       n = n - t2;
     }

     void func_11()
     {
       p = mat3<double>(1,0,0, 0,1,0, 0,0,1);
       cctbx::adptbx::eigensystem<double> n_eigensystem(n);
       vec3<double> nyu = n_eigensystem.values();
       bn = nyu[0];
     }

     void func_12()
     {
       double bm = beta(m, small);
       t1 = sym_mat3<double>(std::min(bm, bn),0,0,0,0,0);
       t_ = t3 + sym_mat3<double>(
          g * (t2 + sym_mat3<double>(p * t1 * p.transpose(), small)) *
          g.transpose(), small);
       MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(t_,       small));
       bool eval1 = cctbx::adptbx::is_positive_definite(m_cp-t_, small3);
       bool eval2 = cctbx::adptbx::is_positive_definite(n_cp-t_, small3);
       if(!eval1 || !eval2) show_all();
       MMTBX_ASSERT(eval1);
       MMTBX_ASSERT(eval2);
     }

     double alpha(sym_mat3<double> A)
      {
        double a11 = A[0];
        double a22 = A[1];
        double a33 = A[2];
        double a12 = A[3];
        double a13 = A[4];
        double a23 = A[5];
        double a21 = a12;
        double a31 = a13;
        double a32 = a23;
        double A11 = a22*a33 - a23*a32;
        double A22 = a11*a33 - a13*a31;
        double A33 = a11*a22 - a12*a21;
        double det_A = A.determinant();
        double arg = (A11 + A22) * (A11 + A22) - 4. * a33 * det_A;
        if(std::abs(arg)<small) arg = 0.0;
        MMTBX_ASSERT(a33 != 0.0);
        return (A11 + A22 - std::sqrt(arg)) / (2. * a33);
      }

     double beta(sym_mat3<double> A, double tol)
     {
       double result = 0.0;
       double a11   = zero_scalar(A[0]);
       double a22   = zero_scalar(A[1]);
       double a33   = zero_scalar(A[2]);
       double a12   = zero_scalar(A[3]);
       double a13   = zero_scalar(A[4]);
       double a23   = zero_scalar(A[5]);
       double a21   = zero_scalar(a12);
       double a31   = zero_scalar(a13);
       double a32   = zero_scalar(a23);
       double A11   = zero_scalar(a22*a33 - a23*a32);
       double A22   = zero_scalar(a11*a33 - a13*a31);
       double A33   = zero_scalar(a11*a22 - a12*a21);
       double det_A = zero_scalar(A.determinant());
       if(det_A > 0 && A11 != 0.0) {
          result = det_A / A11;
       }
       else if(det_A == 0.0 && A11 > 0) {
          result = 0.0;
       }
       else if(det_A == 0.0 && A11 == 0.0 && (A22+A33) >0) {
          result = (A22+A33) / (a22+a33);
       }
       else if(det_A == 0.0 && A11 == 0.0 && a22 == 0.0 && a33 == 0.0) {
          result = a11;
       }
       else if(det_A==0.0 && A11==0.0 && A22==0.0 && A33==0.0 && (a22+a33)>0) {
          result = 0.0;
       }
       // exception case
       else if(det_A == 0.0 || A11 == 0.0) {
          result = 0.0;
       }
       else {
          return 0.0;
          //show_all();
          //std::cout<<"det_A "<<det_A<<std::endl;
          //std::cout<<"A11 "<<A11<<std::endl;
          //std::cout<<"A22 "<<A22<<std::endl;
          //std::cout<<"A33 "<<A33<<std::endl;
          //std::cout<<"A22+A33 "<<A22+A33<<std::endl;
          //std::cout<<"a22 "<<a22<<std::endl;
          //std::cout<<"a33 "<<a33<<std::endl;
          //std::cout<<"a22+a33 "<<a22+a33<<std::endl;
          //throw error("None of possible choices.");
       }
       return result;
     }

     void show_mat(sym_mat3<double> x)
     {
        std::cout<<"matrix= "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<" "<<x[4]<<" "<<x[5]<<std::endl;
     }

     void show_vec(vec3<double> x)
     {
        std::cout<<"vec= "<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
     }

     void show_all() {
       std::cout<<"***********************************************"<<std::endl;
       std::cout<<"Start values: "<<std::endl;
       cctbx::adptbx::eigensystem<double> m_eigensystem_cp(m_cp);
       cctbx::adptbx::eigensystem<double> n_eigensystem_cp(n_cp);
       vec3<double> nyu_cp = n_eigensystem_cp.values();
       vec3<double> myu_cp = m_eigensystem_cp.values();
       std::cout<<"M= "<<m_cp[0]<<" "<<m_cp[1]<<" "<<m_cp[2]<<" "<<m_cp[3]<<" "<<m_cp[4]<<" "<<m_cp[5]<<std::endl;
       std::cout<<"N= "<<n_cp[0]<<" "<<n_cp[1]<<" "<<n_cp[2]<<" "<<n_cp[3]<<" "<<n_cp[4]<<" "<<n_cp[5]<<std::endl;
       std::cout<<"myu= "<<myu_cp[0]<<" "<<myu_cp[1]<<" "<<myu_cp[2]<<std::endl;
       std::cout<<"nyu= "<<nyu_cp[0]<<" "<<nyu_cp[1]<<" "<<nyu_cp[2]<<std::endl;
       std::cout<<"Current values: "<<std::endl;
       cctbx::adptbx::eigensystem<double> m_eigensystem(m);
       cctbx::adptbx::eigensystem<double> n_eigensystem(n);
       vec3<double> nyu = n_eigensystem.values();
       vec3<double> myu = m_eigensystem.values();
       std::cout<<"myu= "<<myu[0]<<" "<<myu[1]<<" "<<myu[2]<<std::endl;
       std::cout<<"nyu= "<<nyu[0]<<" "<<nyu[1]<<" "<<nyu[2]<<std::endl;
       std::cout<<"M= "<<m[0]<<" "<<m[1]<<" "<<m[2]<<" "<<m[3]<<" "<<m[4]<<" "<<m[5]<<std::endl;
       std::cout<<"N= "<<n[0]<<" "<<n[1]<<" "<<n[2]<<" "<<n[3]<<" "<<n[4]<<" "<<n[5]<<std::endl;
       std::cout<<"T= "<<t_[0]<<" "<<t_[1]<<" "<<t_[2]<<" "<<t_[3]<<" "<<t_[4]<<" "<<t_[5]<<std::endl;
       std::cout<<"T1= "<<t1[0]<<" "<<t1[1]<<" "<<t1[2]<<" "<<t1[3]<<" "<<t1[4]<<" "<<t1[5]<<std::endl;
       std::cout<<"T2= "<<t2[0]<<" "<<t2[1]<<" "<<t2[2]<<" "<<t2[3]<<" "<<t2[4]<<" "<<t2[5]<<std::endl;
       std::cout<<"T3= "<<t3[0]<<" "<<t3[1]<<" "<<t3[2]<<" "<<t3[3]<<" "<<t3[4]<<" "<<t3[5]<<std::endl;

       cctbx::adptbx::eigensystem<double> m_eigensystem_(m_cp-t_);
       cctbx::adptbx::eigensystem<double> n_eigensystem_(n_cp-t_);
       vec3<double> nyu_ = n_eigensystem_.values();
       vec3<double> myu_ = m_eigensystem_.values();
       std::cout<<"myu(M-T)= "<<myu_[0]<<" "<<myu_[1]<<" "<<myu_[2]<<std::endl;
       std::cout<<"nyu(N-T)= "<<nyu_[0]<<" "<<nyu_[1]<<" "<<nyu_[2]<<std::endl;
       std::cout<<"PD(M)= "<<cctbx::adptbx::is_positive_definite(m, small)<<std::endl;
       std::cout<<"PD(N)= "<<cctbx::adptbx::is_positive_definite(n, small)<<std::endl;


       std::cout<<"branch_0       = "<<branch_0      <<std::endl;
       std::cout<<"branch_1       = "<<branch_1      <<std::endl;
       std::cout<<"branch_1_1     = "<<branch_1_1    <<std::endl;
       std::cout<<"branch_1_2     = "<<branch_1_2    <<std::endl;
       std::cout<<"branch_1_2_1   = "<<branch_1_2_1  <<std::endl;
       std::cout<<"branch_1_2_2   = "<<branch_1_2_2  <<std::endl;
       std::cout<<"branch_1_2_3   = "<<branch_1_2_3  <<std::endl;
       std::cout<<"branch_1_2_3_1 = "<<branch_1_2_3_1<<std::endl;
       std::cout<<"branch_1_2_3_2 = "<<branch_1_2_3_2<<std::endl;
       std::cout<<"***********************************************"<<std::endl;
     }

     sym_mat3<double> zero_sym_mat(sym_mat3<double> x)
     {
       sym_mat3<double> y = x;
       for(int i=0; i<6; i++) {
           if(std::abs(x[i]) < small) y[i] = 0.0;
       }
       return y;
     }

     double zero_scalar(double x)
     {
       if(std::abs(x) < small) return 0.0;
       else return x;
     }

    sym_mat3<double> t() const { return t_; }
    bool branch_0,branch_1,branch_1_1,branch_1_2,branch_1_2_1,branch_1_2_2;
    bool branch_1_2_3,branch_1_2_3_1,branch_1_2_3_2;
    bool get_branch_0()       const {return branch_0      ; }
    bool get_branch_1()       const {return branch_1      ; }
    bool get_branch_1_1()     const {return branch_1_1    ; }
    bool get_branch_1_2()     const {return branch_1_2    ; }
    bool get_branch_1_2_1()   const {return branch_1_2_1  ; }
    bool get_branch_1_2_2()   const {return branch_1_2_2  ; }
    bool get_branch_1_2_3()   const {return branch_1_2_3  ; }
    bool get_branch_1_2_3_1() const {return branch_1_2_3_1; }
    bool get_branch_1_2_3_2() const {return branch_1_2_3_2; }
  protected:
    sym_mat3<double> t1,t2,t3,m,n,m_cp,n_cp;
    mat3<double> p, g;
    double small,small6,small5,small4,small3,a,bn;
};


af::shared<sym_mat3<double> > uaniso_from_tls_one_group(
                                   tlso<double> tls_params,
                                   af::shared<vec3<double> > const& sites_cart,
                                   bool zeroize_trace)
{
    af::shared<sym_mat3<double> > uanisos(sites_cart.size());
    for (std::size_t i=0; i < sites_cart.size(); i++) {
         vec3<double> const& site = sites_cart[i];
         uaniso_from_tls manager(tls_params.t,
                                 tls_params.l,
                                 tls_params.s,
                                 tls_params.origin,
                                 site,
                                 zeroize_trace);
         uanisos[i] = manager.u();
    }
    return uanisos;
}

sym_mat3<double> t_from_u_cart(af::shared<sym_mat3<double> > const& u_cart,
                               double small)
{
    MMTBX_ASSERT(u_cart.size() >= 2);
    int i_seq = -1;
    double ev_min = 1.e+20;
    for(std::size_t i=0; i < u_cart.size(); i++) {
        MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(u_cart[i], 1.e-6));
        cctbx::adptbx::eigensystem<double> eigensystem(u_cart[i]);
        double try_ev = eigensystem.values().min();
        if(try_ev < ev_min) {
           ev_min = try_ev;
           i_seq = i;
        }
    }
    cctbx::adptbx::eigensystem<double> eigensystem_(u_cart[i_seq]);
    MMTBX_ASSERT( eigensystem_.values().min() == ev_min );

    sym_mat3<double> t = common(u_cart[0], u_cart[i_seq]).t();
    if(u_cart.size() > 2) {
       for (std::size_t i=0; i < u_cart.size(); i++) {
         MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(t,         small));
         MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(u_cart[i], small));
         t = common(t, u_cart[i], small).t();
         MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(t,          small));
         MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(u_cart[i]-t,small));
         MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(u_cart[i],  small));
         //for (std::size_t j=0; j < u_cart.size(); j++) {
         //  t = common(t, u_cart[j], small).t();
         //  MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(u_cart[j]-t, small));
         //}
       }
       for (std::size_t i=0; i < u_cart.size(); i++) {
         MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(u_cart[i]-t, small));
       }
    }

    return t;
}

sym_mat3<double> t_from_u_cart(af::shared<double> const& u_iso,
                               double small)
{
    af::shared<sym_mat3<double> > u_cart(u_iso.size());
    for(std::size_t i=0; i < u_cart.size(); i++) {
      u_cart[i] = sym_mat3<double>(u_iso[i], u_iso[i], u_iso[i], 0, 0, 0);
    }

    MMTBX_ASSERT(u_cart.size() >= 2);
    int i_seq = -1;
    double ev_min = 1.e+20;
    for(std::size_t i=0; i < u_cart.size(); i++) {
        MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(u_cart[i], small));
        cctbx::adptbx::eigensystem<double> eigensystem(u_cart[i]);
        double try_ev = eigensystem.values().min();
        if(try_ev < ev_min) {
           ev_min = try_ev;
           i_seq = i;
        }
    }
    cctbx::adptbx::eigensystem<double> eigensystem_(u_cart[i_seq]);
    MMTBX_ASSERT( eigensystem_.values().min() == ev_min );

    sym_mat3<double> t = common(u_cart[0], u_cart[i_seq]).t();
    if(u_cart.size() > 2) {
       for (std::size_t i=0; i < u_cart.size(); i++) {
         MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(t,         small));
         MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(u_cart[i], small));
         t = common(t, u_cart[i], small).t();
         MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(t,          small));
         MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(u_cart[i]-t,small));
         MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(u_cart[i],  small));
         //for (std::size_t j=0; j < u_cart.size(); j++) {
         //  t = common(t, u_cart[j], small).t();
         //  MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(u_cart[j]-t, small));
         //}
       }
       for (std::size_t i=0; i < u_cart.size(); i++) {
         MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(u_cart[i]-t, small));
       }
    }

    return t;
}

class tls_ls_derivative_coefficients {
public:
  af::versa<double, af::flex_grid<> > a;
  af::shared<double> b;
  af::versa<double, af::flex_grid<> > tmp;

  tls_ls_derivative_coefficients(
    vec3<double> const& origin,
    af::shared<vec3<double> > const& sites_cart,
    af::shared<double> const& u_iso)
  :
  a(af::flex_grid<>(10, 10), 0),
  b(10, 0), tmp(af::flex_grid<>(10, 10), 0)
  {
    MMTBX_ASSERT(sites_cart.size() == u_iso.size());
    MMTBX_ASSERT(sites_cart.size() > 0);
    double deg2rad = scitbx::deg_as_rad(1.0);
    double deg2radsq = deg2rad * deg2rad;
    for (std::size_t i=0; i < sites_cart.size(); i++) {
      vec3<double> r = sites_cart[i] - origin;
      double x = r[0];
      double y = r[1];
      double z = r[2];
      //
      double X = (y*y+z*z)/3  * deg2radsq;
      double Y = (x*x+z*z)/3  * deg2radsq;
      double Z = (x*x+y*y)/3  * deg2radsq;
      double W = -2*x*y/3     * deg2radsq;
      double V = -2*x*z/3     * deg2radsq;
      double T = -2*y*z/3     * deg2radsq;
      double S = 2*z/3        * deg2rad;
      double R = 2*y/3        * deg2rad;
      double Q = 2*x/3        * deg2rad;
      //
      const double cv_values[] = {1,X,Y,Z,W,V,T,S,R,Q};
      af::shared<double> cv(cv_values, cv_values+10);
      scitbx::matrix::outer_product(tmp.begin(),cv.const_ref(),cv.const_ref());

      a += tmp;
      b = b + u_iso[i]*cv;
    }
  }
};

double ls_target_from_iso_tls(double const& t,
                              sym_mat3<double> const& l_deg,
                              vec3<double> const& s_deg,
                              vec3<double> const& origin,
                              af::shared<vec3<double> > const& sites,
                              af::shared<double> const& uisos)
  {
    double tg = 0.0;
    double scale_l_and_s = scitbx::deg_as_rad(1.0);
    for(std::size_t i=0; i < sites.size(); i++) {
      vec3<double> const& site = sites[i];
      double diff = uiso_from_tls(t,l_deg,s_deg,origin,site) - uisos[i];
      double two_diff = 2*diff;
      tg += diff*diff;
    }
    return tg;
}

}} // namespace mmtbx::tls

#endif // MMTBX_TLS_H
