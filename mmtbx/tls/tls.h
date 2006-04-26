#ifndef MMTBX_TLS_H
#define MMTBX_TLS_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <vector>
#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>
#include <cmath>
#include <cctbx/adptbx.h>
#include <cctbx/xray/scattering_type_registry.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <mmtbx/error.h>
#include <cctbx/xray/targets.h>

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

class uaniso_from_tls {
public:
  uaniso_from_tls(sym_mat3<double> const& T,
                  sym_mat3<double> const& L_deg,
                  mat3<double> const& S_deg,
                  vec3<double> const& origin,
                  vec3<double> const& site)
  {
   double deg2rad = scitbx::deg_as_rad(1.0);
   double deg2radsq = deg2rad * deg2rad;
   L = L_deg * deg2radsq;
   S = S_deg * deg2rad;
   S[8] = -(S[0]+S[4]); // condition on trace(S) = 0.0
   vec3<double> r = site - origin;
   x = r[0];
   y = r[1];
   z = r[2];
   //xx = x*x;
   //yy = y*y;
   //zz = z*z;
   //xy = x*y;
   //yz = y*z;
   //xz = x*z;
   sym_mat3<double> ALA = L.antisymmetric_tensor_transform(z,-y,x);
   sym_mat3<double> ASSA = sym_mat3<double>(2.*S[3]*z - 2.*S[6]*y,
                                            2.*S[7]*x - 2.*S[1]*z,
                                            2.*S[2]*y - 2.*S[5]*x,
                                            S[6]*x - S[7]*y + (S[4] - S[0])*z,
                                            S[5]*z - S[3]*x + (S[0] - S[8])*y,
                                            S[1]*y - S[2]*z + (S[8] - S[4])*x);
   uaniso = sym_mat3<double>(T + ALA + ASSA);
   //double u11=T[0]+zz*L[1]+yy*L[2]-2.*(yz*L[5]+y*S[6]-z*S[3]);
   //double u12=T[3]-xy*L[2]+xz*L[5]+yz*L[4]-zz*L[3]-z*(S[0]-S[4])+x*S[6]-y*S[7];
   //double u13=T[4]-S[3]*x+y*(S[0]-S[8])+S[5]*z+L[5]*xy-L[1]*xz-L[4]*yy+L[3]*yz;
   //double u22=T[1]+2.*(S[7]*x-S[1]*z-L[4]*xz)+L[0]*zz+L[2]*xx;
   //double u23=T[5]-x*(S[4]-S[8])+S[1]*y-S[2]*z+L[4]*xy-L[0]*yz-L[5]*xx+L[3]*xz;
   //double u33=T[2]-2.*(S[5]*x-S[2]*y)-L[3]*xy+L[0]*yy+L[1]*xx-L[3]*xy;
   //uaniso = sym_mat3<double>(u11,u22,u33,u12,u13,u23);
  }
  double x,y,z,xx,yy,zz,xy,yz,xz;
  sym_mat3<double> u() { return uaniso; }
  sym_mat3<double> Lrad() { return L; }
  mat3<double> Srad() { return S; }
private:
  sym_mat3<double> uaniso;
  sym_mat3<double> L;
  mat3<double> S;
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
         uaniso_from_tls manager(T,L_deg,S_deg,origin,site);
         sym_mat3<double> const& utls = manager.u();
         sym_mat3<double> diff = utls - uanisos[i];
         for (std::size_t k=0; k < diff.size(); k++) {
              tg += diff[k] * diff[k];
         }
         diffs.push_back(diff*2.);
    }
    d_target_d_tls  d_target_d_tls_manager(sites, origin, diffs, false, false);
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


af::shared<sym_mat3<double> > uaniso_from_tls_one_group(
                                   tlso<double> tls_params,
                                   af::shared<vec3<double> > const& sites_cart)
{
    af::shared<sym_mat3<double> > uanisos(sites_cart.size());
    for (std::size_t i=0; i < sites_cart.size(); i++) {
         vec3<double> const& site = sites_cart[i];
         uaniso_from_tls manager(tls_params.t,
                                 tls_params.l,
                                 tls_params.s,
                                 tls_params.origin,
                                 site);
         uanisos[i] = manager.u();
    }
    return uanisos;
}



}} // namespace mmtbx::tls

#endif // MMTBX_TLS_H
