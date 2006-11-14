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
                  vec3<double> const& site_cart)
  {
   double deg2rad = scitbx::deg_as_rad(1.0);
   double deg2radsq = deg2rad * deg2rad;
   L = L_deg * deg2radsq;
   S = S_deg * deg2rad;
   S[8] = -(S[0]+S[4]); // condition on trace(S) = 0.0
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
         uaniso_from_tls manager(T,L_deg,S_deg,origin,site);
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
                                site);
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
                                site);
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
      qg(-1,-1,-1,-1,-1,-1,-1,-1,-1),
      qp(-1,-1,-1,-1,-1,-1,-1,-1,-1),
      q(-1,-1,-1,-1,-1,-1,-1,-1,-1),
      g1(-1,-1,-1), g2(-1,-1,-1), g3(-1,-1,-1),
      myu(0,0,0),
      nyu(0,0,0),
      m(m_),m_cp(m_),
      n(n_),n_cp(n_),
      small(small),small6(1.e-6),
      branch_0(false),branch_1(false),branch_2(false),branch_3(false),
      branch_1_1(false),branch_1_2(false),branch_3_1(false),branch_3_2(false),
      branch_3_3_1(false),branch_3_3_2(false),branch_2_1(false),
      branch_2_2(false),branch_2_3(false),branch_3_3_3(false)
    {
      check_and_flip();
      if(!(myu[2] < nyu[0])) {
         t_ = n;
         branch_0 = true;
      }
      else {
         func_2();
         branch_123();
         if(branch_1) {
            func_1_1();
            func_1_2();
            if(std::abs(nyu[1]) < small) {
               func_1_3();
               branch_1_1 = true;
            }
            else {
               branch_1_2 = true;
               goto jump;
            }
         }
         if(branch_2) {
jump:
            func_2_1();
            func_2_2();
         }
         if(branch_3) {
            func_3_1();
            if(std::abs(nyu[1]) < small) {
               branch_3_1 = true;
               func_3_3();
            }
            else {
               branch_3_2 = true;
               func_3_2();
               func_3_3();
            }
         }
      }
      bool zz1 = cctbx::adptbx::is_positive_definite(t_,      small6*10.);
      bool zz2 = cctbx::adptbx::is_positive_definite(m_cp-t_, small6*10.);
      bool zz3 = cctbx::adptbx::is_positive_definite(n_cp-t_, small6*10.);
      if(!zz1 || !zz2 || !zz3) show_all();
      MMTBX_ASSERT( cctbx::adptbx::is_positive_definite(t_,      small6*10.));//XXX
      MMTBX_ASSERT( cctbx::adptbx::is_positive_definite(m_cp-t_, small6*10.));//XXX
      MMTBX_ASSERT( cctbx::adptbx::is_positive_definite(n_cp-t_, small6*10.));//XXX
    }

    void check_and_flip()
     {
       cctbx::adptbx::eigensystem<double> m_eigensystem(m);
       cctbx::adptbx::eigensystem<double> n_eigensystem(n);
       myu = m_eigensystem.values();
       nyu = n_eigensystem.values();
       MMTBX_ASSERT( cctbx::adptbx::is_positive_definite(myu, small6) );
       MMTBX_ASSERT( cctbx::adptbx::is_positive_definite(nyu, small6) );
       MMTBX_ASSERT( myu[0]>=myu[1] && myu[1]>=myu[2] && myu[2]>0. );
       MMTBX_ASSERT( nyu[0]>=nyu[1] && nyu[1]>=nyu[2] && nyu[2]>0. );
       if((myu[2]<nyu[2])||(std::abs(myu[2]-nyu[2])<small && myu[1]<nyu[1])||
          (std::abs(myu[2]-nyu[2])<small&&std::abs(myu[1]-nyu[1])<small &&
          myu[0]<nyu[0])) {
          sym_mat3<double> tmp = m;
          m = n;
          n = tmp;
          vec3<double> tmp_ = myu;
          myu = nyu;
          nyu = tmp_;
       }
     }

     void func_2()
      {
        t3 = sym_mat3<double>(nyu[2],nyu[2],nyu[2],0,0,0);
        m = m - t3;
        n = n - t3;
        myu = vec3<double>(myu[0]-nyu[2], myu[1]-nyu[2], myu[2]-nyu[2]);
        nyu = vec3<double>(nyu[0]-nyu[2], nyu[1]-nyu[2], 0);
        cctbx::adptbx::eigensystem<double> n_eigensystem(n);
        MMTBX_ASSERT(n_eigensystem.values().const_ref().all_approx_equal(
                                                      nyu.const_ref(), small));
        g1 = n_eigensystem.vectors(0);
        g2 = n_eigensystem.vectors(1);
        g3 = n_eigensystem.vectors(2);
      }

     void branch_123()
      {
        if(!(std::abs(myu[2]) < small)) branch_1 = true;
        else {
           double if_zero = g3 * (m * g3);
           if(std::abs(if_zero) < small) branch_3 = true;
           else branch_2 = true;
        }
        MMTBX_ASSERT(branch_1 || branch_2 || branch_3);
      }

     void func_1_1()
      {
        n = sym_mat3<double>(nyu[0],nyu[1],0,0,0,0);
        double g1_n = g1.length();
        double g2_n = g2.length();
        double g3_n = g3.length();
        qg = mat3<double>(g1[0]/g1_n,g2[0]/g2_n,g3[0]/g3_n,
                          g1[1]/g1_n,g2[1]/g2_n,g3[1]/g3_n,
                          g1[2]/g1_n,g2[2]/g2_n,g3[2]/g3_n);
        m = sym_mat3<double>(qg.inverse() * m * qg, small);
      }

     void func_1_2()
      {
        double det_m = m.determinant();
        double det_m11 = m[1]*m[2] - m[5]*m[5];
        double det_m22 = m[0]*m[2] - m[4]*m[4];
        double det_m33 = m[0]*m[1] - m[3]*m[3];
        double arg = (det_m11+det_m22)*(det_m11+det_m22) - 4.*m[2]* det_m;
        if(arg < 0.0) arg = 0.0;
        MMTBX_ASSERT(std::abs(m[2]) >= small);
        double alpha_1 = (det_m11 + det_m22 - std::sqrt(arg)) / (2.*m[2]);
        MMTBX_ASSERT(alpha_1 > 0.0);
        double alpha = std::min(alpha_1, nyu[1]);
        t2 = sym_mat3<double>(alpha,alpha,0,0,0,0);
        m = m - t2;
        n = n - t2;
        nyu = vec3<double>(nyu[0]-alpha, nyu[1]-alpha, 0);
        cctbx::adptbx::eigensystem<double> m_eigensystem(m);
        cctbx::adptbx::eigensystem<double> n_eigensystem(n);
        MMTBX_ASSERT(n_eigensystem.values().const_ref().all_approx_equal(
                                                      nyu.const_ref(), small));
        myu = m_eigensystem.values();
      }

     void func_1_3()
      {
        MMTBX_ASSERT(myu[0]>=myu[1]&&myu[1]>=myu[2]&&
                    (std::abs(myu[2])<small || myu[2] > 0.0));
        MMTBX_ASSERT(nyu[0]>=nyu[1]&&std::abs(nyu[1]-nyu[2])<small&&
                     std::abs(nyu[2])<small);
        double det_m = m.determinant();
        double det_m11 = m[1]*m[2] - m[5]*m[5];
        double det_m22 = m[0]*m[2] - m[4]*m[4];
        double det_m33 = m[0]*m[1] - m[3]*m[3];
        double beta_1 = get_beta_1(det_m, det_m11, det_m22, det_m33, m[0],
                                                         m[1], m[2], small);
        double beta = std::min(beta_1, nyu[0]);
        t1 = sym_mat3<double>(beta,0,0,0,0,0);
        t_ = t3 + sym_mat3<double>(qg * (t2 + t1) * qg.inverse(), small);
        MMTBX_ASSERT( cctbx::adptbx::is_positive_definite(t_,      small6));
        MMTBX_ASSERT( cctbx::adptbx::is_positive_definite(m_cp-t_, small6));
        MMTBX_ASSERT( cctbx::adptbx::is_positive_definite(n_cp-t_, small6));
      }

     void func_2_1()
      {
        if(branch_2) {
           MMTBX_ASSERT(myu[0]>=myu[1]&&myu[1]>=myu[2]&&std::abs(myu[2])<small);
           MMTBX_ASSERT(nyu[0]>=nyu[1]&&nyu[1]>=nyu[2]&&std::abs(nyu[2])<small);
        }
        else {
           MMTBX_ASSERT(myu[0]>=myu[1]&&myu[1]>=myu[2]&&std::abs(myu[2])<small);
           MMTBX_ASSERT(nyu[0]>=nyu[1]&&nyu[1]>=nyu[2]&&
                        std::abs(nyu[1])>small&&std::abs(nyu[2])<small);
        }
        cctbx::adptbx::eigensystem<double> n_eigensystem(n);
        cctbx::adptbx::eigensystem<double> m_eigensystem(m);
        MMTBX_ASSERT(m_eigensystem.values().const_ref().all_approx_equal(
                                                      myu.const_ref(), small));
        vec3<double> e1 = m_eigensystem.vectors(0);
        vec3<double> e2 = m_eigensystem.vectors(1);
        vec3<double> e3 = m_eigensystem.vectors(2);
        g1 = n_eigensystem.vectors(0);
        g2 = n_eigensystem.vectors(1);
        g3 = n_eigensystem.vectors(2);
        MMTBX_ASSERT((e3-g3).length() > small6);
        double lambda = - (e3 * g3);
        vec3<double> p2 = e3 + g3*lambda;
        vec3<double> p3 = g3;
        p2 = p2 / p2.length();
        p3 = p3 / p3.length();
        vec3<double> p1 = p2.cross(p3);
        p1 = p1 / p1.length();
        qp = mat3<double>(p1[0],p2[0],p3[0],
                          p1[1],p2[1],p3[1],
                          p1[2],p2[2],p3[2]);
        //qp = zero_mat(qp);
        m = sym_mat3<double>(qp.inverse() * m * qp, small);
        n = sym_mat3<double>(qp.inverse() * n * qp, small);
        //zero_all();
      }

     void func_2_2()
      {
        if(branch_2) {
           MMTBX_ASSERT(std::abs(n[2])<small6 && std::abs(n[4])<small6 &&
                        std::abs(n[5])<small6);
           MMTBX_ASSERT(!branch_1_2 && !branch_1_1 && !branch_1);
           n = sym_mat3<double>(n[0],n[1],0,n[3],0,0);
        }
        cctbx::adptbx::eigensystem<double> n_eigensystem(n);
        nyu = n_eigensystem.values();
        //zero_all();
        double det_m   = m.determinant();
        double det_m11 = m[1]*m[2] - m[5]*m[5];
        double det_m22 = m[0]*m[2] - m[4]*m[4];
        double det_m33 = m[0]*m[1] - m[3]*m[3];
        double beta_1 = get_beta_1(det_m, det_m11, det_m22, det_m33, m[0],
                                                            m[1], m[2], small);
        double gamma_2 = 999999.0;
        double det_n33 = n[0]*n[1] - n[3]*n[3];
        if(det_n33 > small) {
           branch_2_1 = true;
           gamma_2 = det_n33 / n[1];
        }
        else {
           MMTBX_ASSERT(std::abs(det_n33) < small);
           if(n[1] > small) {
              branch_2_2 = true;
              gamma_2 = 0.0;
           }
           else {
              branch_2_3 = true;
              MMTBX_ASSERT(std::abs(n[1]) < small);
              gamma_2 = n[0];
           }
        }
        MMTBX_ASSERT(gamma_2 != 999999.0);
        double tau = std::min(beta_1, gamma_2);
        t1 = sym_mat3<double>(tau,0,0,0,0,0);
        if(branch_2) {
           //zero_all();
           t_ = t3 + sym_mat3<double>(qp * t1 * qp.inverse(), small);
           MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(t_,      small6));
           MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(m_cp-t_, 1.e-4));
           MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(n_cp-t_, small6));
        }
        else {
           MMTBX_ASSERT(branch_1_2);
           //zero_all();
           t_ = t3 + sym_mat3<double>(qg * (t2 + sym_mat3<double>(
                          qp * t1 * qp.inverse(),small)) * qg.inverse(),small);
           MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(t_,      small6));
           MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(m_cp-t_, small6));
           MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(n_cp-t_, small6));
        }
      }

     void func_3_1()
      {
        MMTBX_ASSERT(myu[0]>=myu[1]&&myu[1]>=myu[2]&&std::abs(myu[2])<small);
        MMTBX_ASSERT(nyu[0]>=nyu[1]&&nyu[1]>=nyu[2]&&std::abs(nyu[2])<small);
        cctbx::adptbx::eigensystem<double> n_eigensystem(n);
        vec3<double> e1 = n_eigensystem.vectors(0);
        vec3<double> e2 = n_eigensystem.vectors(1);
        vec3<double> e3 = n_eigensystem.vectors(2);
        MMTBX_ASSERT(g3.const_ref().all_approx_equal(
                     e3.const_ref(),small));
        MMTBX_ASSERT(nyu.const_ref().all_approx_equal(
                     n_eigensystem.values().const_ref(),small));
        //zero_all();
        n = sym_mat3<double>(nyu[0],nyu[1],0,0,0,0);
        double e1_n = e1.length();
        double e2_n = e2.length();
        double e3_n = e3.length();
        q = mat3<double>(e1[0]/e1_n,e2[0]/e2_n,e3[0]/e3_n,
                         e1[1]/e1_n,e2[1]/e2_n,e3[1]/e3_n,
                         e1[2]/e1_n,e2[2]/e2_n,e3[2]/e3_n);
        //q = zero_mat(q);
        m = sym_mat3<double>(q.inverse() * m * q, small);
        if(!(std::abs(m[4])<1.e-4 && std::abs(m[4])<1.e-4 &&
                                                       std::abs(m[5])<1.e-4)) {
           show_all();
        }
        //MMTBX_ASSERT(std::abs(m[4])<1.e-4 && std::abs(m[4])<1.e-4 &&
        //             std::abs(m[5])<1.e-4);
        m = sym_mat3<double>(m[0],m[1],0,m[3],0,0);
      }

     void func_3_2()
      {
        //MMTBX_ASSERT(nyu[0]>=nyu[1]&&std::abs(nyu[1])>small&&myu[1]>=nyu[1]);
        MMTBX_ASSERT(nyu[0]>=nyu[1]);
        MMTBX_ASSERT(std::abs(nyu[1])>small);
        //MMTBX_ASSERT(myu[1]>=nyu[1]); //XXX
        MMTBX_ASSERT((std::abs(myu[1]-nyu[1]) < small6) || (myu[1]>=nyu[1]));
        cctbx::adptbx::eigensystem<double> m_eigensystem(m);
        myu = m_eigensystem.values();
        t2 = sym_mat3<double>(nyu[1],nyu[1],0,0,0,0);
        m = m - t2;
        n = n - t2;
        myu = vec3<double>(myu[0]-nyu[1], myu[1]-nyu[1], 0);
        nyu = vec3<double>(nyu[0]-nyu[1], 0, 0);
      }

     void func_3_3()
      {
        double gamma_2 = 0.0;
        double det_m33 = m[0]*m[1] - m[3]*m[3];
        if(det_m33 > small) {
           branch_3_3_1 = true;
           gamma_2 = det_m33 / m[1];
        }
        else {
           MMTBX_ASSERT(std::abs(det_m33) < small);
           if(m[1] > small) {
              branch_3_3_2 = true;
              gamma_2 = 0.0;
           }
           else {
              branch_3_3_3 = true;
              MMTBX_ASSERT(std::abs(m[1]) < small);
              gamma_2 = m[0];
           }
        }


        //if(std::abs(myu[1]) < small) {
        //   branch_3_3_1 = true;
        //   MMTBX_ASSERT(myu[0]>=myu[1]&&std::abs(myu[1])<small);
        //   gamma_2 = m[0];
        //}
        //else {
        //   branch_3_3_2 = true;
        //   MMTBX_ASSERT(myu[0]>=myu[1]&&std::abs(myu[1])>small);
        //   double det_m33 = m[0]*m[1] - m[3]*m[3];
        //   gamma_2 = det_m33 / m[1];
        //}


        MMTBX_ASSERT(nyu[0]>=nyu[1]&&std::abs(nyu[1])<small);
        double gamma = std::min(gamma_2, nyu[0]);
        t1 = sym_mat3<double>(gamma,0,0,0,0,0);
        //zero_all();
        if(branch_3_1) {
           MMTBX_ASSERT(t2 == sym_mat3<double>(-1,-1,-1,-1,-1,-1));
           t2 = sym_mat3<double>(0,0,0,0,0,0);
        }
        t_ = t3 + sym_mat3<double>(q * (t2 + t1) * q.inverse(), small);
        MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(t_,      small6));
        MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(m_cp-t_, 1.e-4)); //XXX
        MMTBX_ASSERT(cctbx::adptbx::is_positive_definite(n_cp-t_, 1.e-4)); //XXX
      }

     double get_beta_1(double det_m,
                       double det_m11,
                       double det_m22,
                       double det_m33,
                       double m11,
                       double m22,
                       double m33,
                       double small)
     {
       double result = 2679941.0;
       if(std::abs(det_m11) < small) det_m = 0.0;

       if(det_m > 0.0) {
          result = det_m / det_m11;
       }
       else if (det_m11 > 0.0) {
          result = 0.0;
       }
       else if (m22+m33 > 0.0) {
          result = (det_m22+det_m33)/(m22+m33);
       }
       else {
          result = m11;
       }
       if(result < 0.0) result = 0.0;
       MMTBX_ASSERT(result != 2679941.0);
       return result;
     }

     void show_all() {
       std::cout<<"***"<<std::endl;
       std::cout<<" "<<std::endl;
       std::cout<<"M(start)= "<<m_cp[0]<<" "<<m_cp[1]<<" "<<m_cp[2]<<" "<<
                                m_cp[3]<<" "<<m_cp[4]<<" "<<m_cp[5]<<std::endl;
       std::cout<<"N(start)= "<<n_cp[0]<<" "<<n_cp[1]<<" "<<n_cp[2]<<" "<<
                                n_cp[3]<<" "<<n_cp[4]<<" "<<n_cp[5]<<std::endl;
       cctbx::adptbx::eigensystem<double> m_eigensystem_s(m_cp);
       cctbx::adptbx::eigensystem<double> n_eigensystem_s(n_cp);
       vec3<double> nvs = n_eigensystem_s.values();
       vec3<double> mvs = m_eigensystem_s.values();
       std::cout<<"myu_start= "<<mvs[0]<<" "<<mvs[1]<<" "<<mvs[2]<<std::endl;
       std::cout<<"nyu_start= "<<nvs[0]<<" "<<nvs[1]<<" "<<nvs[2]<<std::endl;
       std::cout<<"T1= "<<t1[0]<<" "<<t1[1]<<" "<<t1[2]<<" "<<t1[3]<<" "<<
                          t1[4]<<" "<<t1[5]<<std::endl;
       std::cout<<"T2= "<<t2[0]<<" "<<t2[1]<<" "<<t2[2]<<" "<<t2[3]<<" "<<
                          t2[4]<<" "<<t2[5]<<std::endl;
       std::cout<<"T3= "<<t3[0]<<" "<<t3[1]<<" "<<t3[2]<<" "<<t3[3]<<" "<<
                          t3[4]<<" "<<t3[5]<<std::endl;
       std::cout<<"T=  "<<t_[0]<<" "<<t_[1]<<" "<<t_[2]<<" "<<t_[3]<<" "<<
                          t_[4]<<" "<<t_[5]<<std::endl;
       std::cout<<"M(current)= "<<m[0]<<" "<<m[1]<<" "<<m[2]<<" "<<m[3]<<" "<<
                                  m[4]<<" "<<m[5]<<std::endl;
       std::cout<<"N(current)= "<<n[0]<<" "<<n[1]<<" "<<n[2]<<" "<<n[3]<<" "<<
                                  n[4]<<" "<<n[5]<<std::endl;
       std::cout<<"myu= "<<myu[0]<<" "<<myu[1]<<" "<<myu[2]<<std::endl;
       std::cout<<"nyu= "<<nyu[0]<<" "<<nyu[1]<<" "<<nyu[2]<<std::endl;
       std::cout<<"branch_0     = "<<branch_0     <<std::endl;
       std::cout<<"branch_1     = "<<branch_1     <<std::endl;
       std::cout<<"branch_2     = "<<branch_2     <<std::endl;
       std::cout<<"branch_3     = "<<branch_3     <<std::endl;
       std::cout<<"branch_1_1   = "<<branch_1_1   <<std::endl;
       std::cout<<"branch_1_2   = "<<branch_1_2   <<std::endl;
       std::cout<<"branch_3_1   = "<<branch_3_1   <<std::endl;
       std::cout<<"branch_3_2   = "<<branch_3_2   <<std::endl;
       std::cout<<"branch_3_3_1 = "<<branch_3_3_1 <<std::endl;
       std::cout<<"branch_3_3_2 = "<<branch_3_3_2 <<std::endl;
       std::cout<<"branch_3_3_3 = "<<branch_3_3_3 <<std::endl;
       std::cout<<"branch_2_1   = "<<branch_2_1   <<std::endl;
       std::cout<<"branch_2_2   = "<<branch_2_2   <<std::endl;
       std::cout<<"branch_2_3   = "<<branch_2_3   <<std::endl;
       cctbx::adptbx::eigensystem<double> m_eigensystem(m_cp-t_);
       cctbx::adptbx::eigensystem<double> n_eigensystem(n_cp-t_);
       vec3<double> nv = n_eigensystem.values();
       vec3<double> mv = m_eigensystem.values();
       std::cout<<"myu-t= "<<mv[0]<<" "<<mv[1]<<" "<<mv[2]<<std::endl;
       std::cout<<"nyu-t= "<<nv[0]<<" "<<nv[1]<<" "<<nv[2]<<std::endl;
       std::cout<<"***"<<std::endl;
     }

     bool all_zero_m3(sym_mat3<double> x)
     {
       return std::abs(x[0])<small&&std::abs(x[1])<small&&
              std::abs(x[2])<small&&std::abs(x[3])<small&&
              std::abs(x[4])<small&&std::abs(x[5])<small;
     }

     mat3<double> zero_mat(mat3<double> x)
     {
       for(int i=0; i<9; i++) {
           if(std::abs(x[i]) < small) x[i] = 0.0;
       }
       return x;
     }

     void zero_all()
     {
       for(int i=0; i<3; i++) {
           if(std::abs(myu[i]) < small) myu[i] = 0.0;
           if(std::abs(nyu[i]) < small) nyu[i] = 0.0;
       }
       for(int i=0; i<6; i++) {
           if(std::abs(m[i]) < small) m[i] = 0.0;
           if(std::abs(n[i]) < small) n[i] = 0.0;
           if(std::abs(t1[i]) < small) t1[i] = 0.0;
           if(std::abs(t2[i]) < small) t2[i] = 0.0;
           if(std::abs(t3[i]) < small) t3[i] = 0.0;
       }
     }

    sym_mat3<double> t() const { return t_; }
    bool branch_0,branch_1,branch_2,branch_3,branch_1_1,branch_1_2;
    bool branch_3_1,branch_3_2,branch_3_3_1,branch_3_3_2,branch_3_3_3;
    bool branch_2_1,branch_2_2,branch_2_3;
    bool get_branch_0()   const {return branch_0  ; }
    bool get_branch_1()   const {return branch_1  ; }
    bool get_branch_2()   const {return branch_2  ; }
    bool get_branch_3()   const {return branch_3  ; }
    bool get_branch_1_1() const {return branch_1_1; }
    bool get_branch_1_2() const {return branch_1_2; }
    bool get_branch_3_1() const {return branch_3_1; }
    bool get_branch_3_2() const {return branch_3_2; }
    bool get_branch_3_3_1() const {return branch_3_3_1; }
    bool get_branch_3_3_2() const {return branch_3_3_2; }
    bool get_branch_3_3_3() const {return branch_3_3_3; }
    bool get_branch_2_1() const {return branch_2_1; }
    bool get_branch_2_2() const {return branch_2_2; }
    bool get_branch_2_3() const {return branch_2_3; }
  protected:
    sym_mat3<double> t1,t2,t3,m,n,m_cp,n_cp;
    mat3<double> qg,qp,q;
    vec3<double> g1,g2,g3,myu,nyu;
    double small,small6;
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

}} // namespace mmtbx::tls

#endif // MMTBX_TLS_H
