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

class grads_u_wrt_tls {
public:
  grads_u_wrt_tls(vec3<double> const& origin,
                  vec3<double> const& site)
  {
    // Order is : sym_mat3<double>(a11,a22,a33,a12,a13,a23);
    //                            ( a0, a1, a2, a3, a4, a5)
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
    gT_ = sym_mat3<sym_mat3<double> >(gT0,gT1,gT2,gT3,gT4,gT5);
    gL0 = sym_mat3<double>(nul,  zz,  yy, nul, nul, -yz);
    gL1 = sym_mat3<double>( zz, nul,  xx, nul, -xz, nul);
    gL2 = sym_mat3<double>( yy,  xx, nul, -xy, nul, nul);
    gL3 = sym_mat3<double>(nul, nul, xy2, -zz,  yz,  xz);
    gL4 = sym_mat3<double>(nul, xz2, nul,  yz, -yy,  xy);
    gL5 = sym_mat3<double>(yz2, nul, nul,  xz,  xy, -xx);
    gL_ = sym_mat3<sym_mat3<double> >(gL0,gL1,gL2,gL3,gL4,gL5);
    gS0 = sym_mat3<double>(nul, nul, nul,  -z,   y, nul);
    gS1 = sym_mat3<double>(nul, -z2, nul, nul, nul,   y);
    gS2 = sym_mat3<double>(nul, nul,  y2, nul, nul,  -z);
    gS3 = sym_mat3<double>( z2, nul, nul, nul,  -x, nul);
    gS4 = sym_mat3<double>(nul, nul, nul,   z, nul,  -x);
    gS5 = sym_mat3<double>(nul, nul, -x2, nul,   z, nul);
    gS6 = sym_mat3<double>(-y2, nul, nul,   x, nul, nul);
    gS7 = sym_mat3<double>(nul,  x2, nul,  -y, nul, nul);
    gS8 = sym_mat3<double>(nul, nul, nul, nul,  -y, x);
    gS_ = mat3<sym_mat3<double> >(gS0,gS1,gS2,gS3,gS4,gS5,gS6,gS7,gS8);
  }
  sym_mat3<sym_mat3<double> > d_u_d_T() { return gT_; }
  sym_mat3<sym_mat3<double> > d_u_d_L() { return gL_; }
  mat3<sym_mat3<double> >     d_u_d_S() { return gS_; }
private:
  double x,y,z,xx,yy,zz,xy,yz,xz,xy2,xz2,yz2,z2,y2,x2;
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
   xx = x*x;
   yy = y*y;
   zz = z*z;
   xy = x*y;
   yz = y*z;
   xz = x*z;
   double u11=T[0]+zz*L[1]+yy*L[2]-2.*(yz*L[5]+y*S[6]-z*S[3]);
   double u12=T[3]-xy*L[2]+xz*L[5]+yz*L[4]-zz*L[3]-z*(S[0]-S[4])+x*S[6]-y*S[7];
   double u13=T[4]-S[3]*x+y*(S[0]-S[8])+S[5]*z+L[5]*xy-L[1]*xz-L[4]*yy+L[3]*yz;
   double u22=T[1]+2.*(S[7]*x-S[1]*z-L[4]*xz)+L[0]*zz+L[2]*xx;
   double u23=T[5]-x*(S[4]-S[8])+S[1]*y-S[2]*z+L[4]*xy-L[0]*yz-L[5]*xx+L[3]*xz;
   double u33=T[2]-2.*(S[5]*x-S[2]*y)-L[3]*xy+L[0]*yy+L[1]*xx-L[3]*xy;
   uaniso = sym_mat3<double>(u11,u22,u33,u12,u13,u23);
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


class tls_from_uaniso_target_and_grads {
public:
  tls_from_uaniso_target_and_grads(sym_mat3<double> const& T,
                                   sym_mat3<double> const& L_deg,
                                   mat3<double> const& S_deg,
                                   vec3<double> const& origin,
                                   af::shared<vec3<double> > const& sites,
                                   af::shared<sym_mat3<double> > const& uanisos)
  {
    tg = 0.0;
    gT.resize(6,0.0);
    gL.resize(6,0.0);
    gS.resize(9,0.0);
    for (std::size_t i=0; i < sites.size(); i++) {
      vec3<double> const& site = sites[i];
      uaniso_from_tls manager(T,L_deg,S_deg,origin,site);
      sym_mat3<double> const& utls = manager.u();
      sym_mat3<double> diff = utls - uanisos[i];
      for (std::size_t k=0; k < diff.size(); k++) {
           tg += diff[k] * diff[k];
      }
      // grads:
      diff *= 2.0;
      grads_u_wrt_tls grads_u_wrt_tls_manager(origin, site);
      sym_mat3<sym_mat3<double> > d_u_d_T = grads_u_wrt_tls_manager.d_u_d_T();
      sym_mat3<sym_mat3<double> > d_u_d_L = grads_u_wrt_tls_manager.d_u_d_L();
      mat3<sym_mat3<double> >     d_u_d_S = grads_u_wrt_tls_manager.d_u_d_S();

      for (std::size_t k=0; k < 6; k++) {
           for (std::size_t m=0; m < 6; m++) {
                gT[k] += diff[m] * d_u_d_T[k][m];
                gL[k] += diff[m] * d_u_d_L[k][m];
          }
      }
      for (std::size_t k=0; k < 9; k++) {
           for (std::size_t m=0; m < 6; m++) {
                gS[k] += diff[m] * d_u_d_S[k][m];
          }
      }

    }
  }
  double target() const { return tg; }
  af::shared<double> grad_T() { return gT; }
  af::shared<double> grad_L() { return gL; }
  af::shared<double> grad_S() { return gS; }
private:
  double tg;
  af::shared<double> gT, gL, gS;
};


af::shared<sym_mat3<double> > uaniso_from_tls_domain(
                       sym_mat3<double> const& T,
                       sym_mat3<double> const& L_deg,
                       mat3<double> const& S_deg,
                       vec3<double> const& origin,
                       af::shared<vec3<double> > const& sites)
{
    af::shared<sym_mat3<double> > uanisos(sites.size());
    for (std::size_t i=0; i < sites.size(); i++) {
      vec3<double> const& site = sites[i];
      uaniso_from_tls manager(T,L_deg,S_deg,origin,site);
      uanisos[i] = manager.u();
    }
    return uanisos;
}

/*--------------TESTS AREA-------------------------------------------------*/
//
// s.f. calculation in P1 space group from anisotropic atoms
//
class gT_gL_gS_ {
public:
  gT_gL_gS_(vec3<double> const& origin,
            cctbx::uctbx::unit_cell const& uc,
            vec3<double> const& site,
            cctbx::miller::index<> const& mi)
  {
    gT.resize(6,0.0);
    gL.resize(6,0.0);
    gS.resize(9,0.0);
    mat3<double> a = uc.fractionalization_matrix();
    double x = site[0] - origin[0];
    double y = site[1] - origin[1];
    double z = site[2] - origin[2];
    double deg2rad = scitbx::deg_as_rad(1.0);
    double deg2rad2 = deg2rad * deg2rad;

    double q0 = a[0]*mi[0] + a[3]*mi[1] + a[6]*mi[2];
    double q1 = a[1]*mi[0] + a[4]*mi[1] + a[7]*mi[2];
    double q2 = a[2]*mi[0] + a[5]*mi[1] + a[8]*mi[2];
    double q3 = q2*y - q1*z;
    double q4 = q2*x - q0*z;
    double q5 = q1*x - q0*y;

    gT[0] = q0 * q0;
    gT[1] = q1 * q1;
    gT[2] = q2 * q2;
    gT[3] = q0 * q1 * 2.0;
    gT[4] = q0 * q2 * 2.0;
    gT[5] = q1 * q2 * 2.0;
    gL[0] = q3 * q3;
    gL[1] = q4 * q4;
    gL[2] = q5 * q5;
    gL[3] =-2.0 * q4 * q3;
    gL[4] = 2.0 * q5 * q3;
    gL[5] =-2.0 * q5 * q4;

    gS[0] =-2.0 * q2 * (q1*x - 2.0*q0*y) - 2.0*q0*q1*z;
    gS[1] = 2.0 * q1 * q3;
    gS[2] = 2.0 * q2 * q3;
    gS[3] =-2.0 * q0 * q4;
    gS[4] =-2.0 * q2 * (2.0*q1*x - q0*y) + 2.0*q0*q1*z;
    gS[5] =-2.0 * q2 * q4;
    gS[6] = 2.0 * q0 * q5;
    gS[7] = 2.0 * q1 * q5;
    gS[8] = 0.0;

    for(std::size_t i=0; i <= 5; i++) {
      gL[i] *= deg2rad2;
    }
    for(std::size_t i=0; i <= 8; i++) {
      gS[i] *= deg2rad;
    }
  }
    af::shared<double> grad_T() { return gT; }
    af::shared<double> grad_L() { return gL; }
    af::shared<double> grad_S() { return gS; }
private:
    af::shared<double> gT, gL, gS;
};




class gT_gL_gS {
public:
  gT_gL_gS(vec3<double> const& origin,
           cctbx::uctbx::unit_cell const& uc,
           vec3<double> const& site,
           cctbx::miller::index<> const& mi)
  {
    gT.resize(6,0.0);
    gL.resize(6,0.0);
    gS.resize(9,0.0);
    mat3<double> a = uc.fractionalization_matrix();
    double x = site[0] - origin[0];
    double y = site[1] - origin[1];
    double z = site[2] - origin[2];
    double deg2rad = scitbx::deg_as_rad(1.0);
    double deg2rad2 = deg2rad * deg2rad;

    double q0 = a[0]*mi[0] + a[3]*mi[1] + a[6]*mi[2];
    double q1 = a[1]*mi[0] + a[4]*mi[1] + a[7]*mi[2];
    double q2 = a[2]*mi[0] + a[5]*mi[1] + a[8]*mi[2];
    double q3 = q2*y - q1*z;
    double q4 = q2*x - q0*z;
    double q5 = q1*x - q0*y;
    gT[0] = q0 * q0;
    gT[1] = q1 * q1;
    gT[2] = q2 * q2;
    gT[3] = q0 * q1 * 2.;
    gT[4] = q0 * q2 * 2.;
    gT[5] = q1 * q2 * 2.;
    gL[0] = q3 * q3;
    gL[1] = q4 * q4;
    gL[2] = q5 * q5;
    gL[3] =-q4 * q3 * 2.;
    gL[4] = q5 * q3 * 2.;
    gL[5] =-q5 * q4 * 2.;
    gS[0] = q0 * q3 * 2.;
    gS[1] = q1 * q3 * 2.;
    gS[2] = q2 * q3 * 2.;
    gS[3] =-q0 * q4 * 2.;
    gS[4] =-q1 * q4 * 2.;
    gS[5] =-q2 * q4 * 2.;
    gS[6] = q0 * q5 * 2.;
    gS[7] = q1 * q5 * 2.;
    gS[8] = q2 * q5 * 2.;

    for(std::size_t i=0; i <= 5; i++) {
      gL[i] *= deg2rad2;
    }
    for(std::size_t i=0; i <= 8; i++) {
      gS[i] *= deg2rad;
    }
  }
    af::shared<double> grad_T() { return gT; }
    af::shared<double> grad_L() { return gL; }
    af::shared<double> grad_S() { return gS; }
private:
    af::shared<double> gT, gL, gS;
};



template <typename FloatType=double>
class params
{
  public:
    scitbx::sym_mat3<FloatType> t;
    scitbx::sym_mat3<FloatType> l;
    scitbx::mat3<FloatType> s;
    params() {
     t.fill(0.0);
     l.fill(0.0);
     s.fill(0.0);
    }
    params(scitbx::sym_mat3<FloatType> const& t_,
           scitbx::sym_mat3<FloatType> const& l_,
           scitbx::mat3<FloatType> const& s_)
    :
      t(t_),l(l_),s(s_)
    {}
};


class tls_xray_target_grads {
public:
  tls_xray_target_grads(af::const_ref<cctbx::miller::index<> > const& hkl,
                        af::const_ref<double> const& fo,
                        af::const_ref< std::complex<double> > const& fc,
                        cctbx::uctbx::unit_cell const& uc,
                        af::const_ref<vec3<double> > const& origin,
                        af::shared<sym_mat3<double> > const& uanisos,
                        cctbx::xray::scattering_type_registry const& scattering_type_registry,
                        af::const_ref<cctbx::xray::scatterer<> > const& scatterers,
                        af::const_ref<int> const& tls_group_flags)
  {
    MMTBX_ASSERT(fo.size() == fc.size() && fo.size() == hkl.size());
    MMTBX_ASSERT(scatterers.size() == uanisos.size());
    int n_groups = af::max(tls_group_flags);
    MMTBX_ASSERT(n_groups == tls_group_flags[tls_group_flags.size()-1]);
    mat3<double> orthogonalization_matrix = uc.orthogonalization_matrix();
    double pi = scitbx::constants::pi;
    double two_pi = pi * 2.0;
    double minus_two_pi = -2.0 * pi * pi;
    af::shared<double> ss = uc.stol_sq(hkl);
    target_ = 0.0;
    gradTLS_.resize(n_groups);

    double num=0.0;
    double denum=0.0;
    double sum_fo_sq=0.0;
    for(std::size_t i=0; i < fo.size(); i++) {
      double fc_abs = std::abs(fc[i]);
      num += fo[i] * fc_abs;
      denum += fc_abs*fc_abs;
      sum_fo_sq += fo[i] * fo[i];
    }
    MMTBX_ASSERT(denum > 0.0);
    MMTBX_ASSERT(sum_fo_sq > 0);
    double scale = num/denum;
    af::shared<std::size_t>
      scattering_type_indices_memory
        = scattering_type_registry.unique_indices(scatterers);
    af::const_ref<std::size_t>
      scattering_type_indices = scattering_type_indices_memory.const_ref();
    for(std::size_t i=0; i < hkl.size(); i++) {
      double delta = fo[i] - scale * std::abs(fc[i]);
      target_ += delta * delta;
      std::complex<double> d_target_d_fc = -2. * scale * delta * std::conj(fc[i]) /
                                                std::abs(fc[i]) / sum_fo_sq;

      cctbx::miller::index<> const& mi = hkl[i];
      af::shared<params<std::complex<double> > > d_fcalc_over_d_tls;
      d_fcalc_over_d_tls.resize(n_groups);
      double tmp1 = 0.0;
      double tmp2 = 0.0;
      double d_star_sq = uc.d_star_sq(mi);
      af::shared<double>
        form_factors_memory
          = scattering_type_registry.unique_form_factors_at_d_star_sq(
              d_star_sq);
      af::const_ref<double> form_factors = form_factors_memory.const_ref();
      for(std::size_t i_seq=0;i_seq<scatterers.size();i_seq++) {
        cctbx::xray::scatterer<> const& scatterer = scatterers[i_seq];
        double f0 = form_factors[scattering_type_indices[i_seq]];
        sym_mat3<double> const& u = scatterer.u_star;
        //sym_mat3<double> const& u = cctbx::adptbx::u_cart_as_u_star(uc, uanisos[i_seq]);
        int tls_group = tls_group_flags[i_seq];
        vec3<double> const& site_frac = scatterer.site;
        vec3<double> const& site_cart = orthogonalization_matrix * site_frac;
        double phase = mi * site_frac * two_pi;
        double dw = cctbx::adptbx::debye_waller_factor_u_star(mi, u);
        double scfact = f0 * dw * minus_two_pi;
        double acalc = scfact*std::cos(phase);
        double bcalc = scfact*std::sin(phase);
        tmp1 += f0*std::cos(phase) * dw;
        tmp2 += f0*std::sin(phase) * dw;
        gT_gL_gS_ manager(origin[tls_group-1], uc, site_cart, mi);
        af::shared<double> gT = manager.grad_T();
        af::shared<double> gL = manager.grad_L();
        af::shared<double> gS = manager.grad_S();
        for (std::size_t p=0; p < 9; p++) {
          if(p < 6) {
            d_fcalc_over_d_tls[tls_group-1].t[p] += std::complex<double> (acalc * gT[p],bcalc * gT[p]);
            d_fcalc_over_d_tls[tls_group-1].l[p] += std::complex<double> (acalc * gL[p],bcalc * gL[p]);
            d_fcalc_over_d_tls[tls_group-1].s[p] += std::complex<double> (acalc * gS[p],bcalc * gS[p]);
          }
          else {
            d_fcalc_over_d_tls[tls_group-1].s[p] += std::complex<double> (acalc * gS[p],bcalc * gS[p]);
          }
        }
      }
      for (std::size_t m=0; m < n_groups; m++) {
        for (std::size_t p=0; p < 9; p++) {
          if(p < 6) {
            gradTLS_[m].t[p] += std::real ( d_target_d_fc * d_fcalc_over_d_tls[m].t[p] );
            gradTLS_[m].l[p] += std::real ( d_target_d_fc * d_fcalc_over_d_tls[m].l[p] );
            gradTLS_[m].s[p] += std::real ( d_target_d_fc * d_fcalc_over_d_tls[m].s[p] );
          }
          else {
            gradTLS_[m].s[p] += std::real ( d_target_d_fc * d_fcalc_over_d_tls[m].s[p] );
          }
        }
      }
    double qq = std::abs(std::abs(std::complex<double>(tmp1,tmp2)) - std::abs(fc[i]));
    if(qq >= 0.0001) cout<<i<<" qq = "<<qq<<endl;
    MMTBX_ASSERT(qq < 0.0001);
    }
    target_ /= sum_fo_sq;
  }

   af::shared<params<double> > gradTLS() { return gradTLS_; };
   double target() { return target_;}
protected:
   af::shared<params<double> > gradTLS_;
   double target_;
};







}} // namespace mmtbx::tls

#endif // MMTBX_TLS_H
