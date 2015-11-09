#ifndef MAX_LIK_H
#define MAX_LIK_H

#include <boost/python/list.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/extract.hpp>

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <vector>
#include <mmtbx/error.h>
#include <cctbx/sgtbx/space_group.h>
#include <scitbx/math/bessel.h>
#include <math.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/shared_reductions.h>

using scitbx::vec3;
namespace af=scitbx::af;
namespace mmtbx {  namespace max_lik {


class wat_dist {
//
// A priori probability distribution of ordered water molecules
// in macromolecular crystals.
// Based on V.Y. Lunin's notes from 10-FEB-2004
// Pavel Afonine / 05-MAR-2004
//

public:
   wat_dist() {}

   void do_wat_dist(double shell,
                af::shared<vec3<double> > const& xyzf,
                af::shared<double> const& atmrad,
                af::shared<std::string> const& element_symbol,
                cctbx::uctbx::unit_cell const& uc,
                cctbx::sgtbx::space_group const& sg,
                vec3<int> const& nxnynz,
                af::shared<int> const& sel_flag,
                double rad,
                int nshells);

   af::versa<double, af::c_grid<3> > data() const { return water_mask_; }

   void as_xplor_map(cctbx::uctbx::unit_cell const& uc,
                     std::string const& outputfile);

   int max_number_of_shells() const { return nshells_; }

private:
   void preparator(cctbx::uctbx::unit_cell const& uc);
   void apply_symmetry_and_make_xyz_within_0_1(
                        af::shared<vec3<double> > const& xyzf,
                        cctbx::sgtbx::space_group const& sg,
                        af::shared<double> const& atmrad,
                        af::shared<int> sel_flag,
                        af::shared<std::string> const& element_symbol);

   af::shared<vec3<double> > xyzf_01;
   af::shared<std::string> element_symbol_;
   af::shared<double> atmrad_01;
   af::shared<int> sel_flag_;
   cctbx::sgtbx::space_group sg_;
   af::versa<double, af::c_grid<3> > water_mask_;
   void set_shells(cctbx::uctbx::unit_cell const& uc, int nshells);
   int NX,NY,NZ,nshells_;
   double rad_;
   double shell,rs[12][240];
   double asq,bsq,csq,tabc,tacc,tbcc,as,bs,cs;
   double xshell,yshell,zshell,xshrink,yshrink,zshrink;
   double accessible_surf_fract,contact_surf_fract;
   inline int nint(double);
};
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class alpha_beta_est {
public:
  alpha_beta_est(boost::python::list const& fo_sets_list,
                 boost::python::list const& fm_sets_list,
                 boost::python::list const& hkl_sets_list,
                 boost::python::list const& epsilons_sets_list,
                 cctbx::sgtbx::space_group const& sg)
     {
       int len_list = boost::python::len(fo_sets_list);
       af::shared<double> A_in_zones = af::shared<double> (len_list);
       af::shared<double> B_in_zones = af::shared<double> (len_list);
       af::shared<double> topt = af::shared<double> (len_list);
       alpha_ = af::shared<double> (len_list);
       beta_  = af::shared<double> (len_list);

       for(std::size_t i=0; i < len_list; i++) {
         boost::python::extract<af::shared<double> > elem_proxy_1(fo_sets_list[i]);
         fo = elem_proxy_1();
         boost::python::extract<af::shared<double> > elem_proxy_2(fm_sets_list[i]);
         fm = elem_proxy_2();
         boost::python::extract<af::shared<cctbx::miller::index<> > > elem_proxy_3(hkl_sets_list[i]);
         hkl = elem_proxy_3();
         boost::python::extract<af::shared<double> > elem_proxy_4(epsilons_sets_list[i]);
         epsilons = elem_proxy_4();
         A_B_topt_est(fo,fm,hkl,epsilons,sg,A_in_zones[i],B_in_zones[i],topt[i]); //going around in a circle for external access
       }
       topt = smooth(topt);
       alpha_beta_in_zones(A_in_zones, B_in_zones, topt);
     }

    void  A_B_topt_est(af::shared<double> const& fo,
                 af::shared<double> const& fm,
                 af::shared<cctbx::miller::index<> > const& hkl,
                 af::shared<double> const& epsilons,
                 cctbx::sgtbx::space_group const& sg,
                 double& ext_A,
                 double& ext_B,
                 double& ext_topt
                 )
     {
         MMTBX_ASSERT(fo.size() > 0 && fm.size() > 0);
         MMTBX_ASSERT(fo.size() == fm.size());
         MMTBX_ASSERT(fo.size() == hkl.size());
         eps = epsilons;
         cf = sg.is_centric(hkl.const_ref());
         A_B_C_D_omega();
         double topt;
         if      (OMEGAi <= 0.0)       topt = 0.0;
         else if (wi/(A*B) <= 3.0E-07) topt = 1.0e+10;
         else {
                                       topt = solvm();
         }
         ext_A = A;
         ext_B = B;
         ext_topt = topt;
     }

  void A_B_C_D_omega()
     {
       SUMwj = 0.0, A = 0.0, B = 0.0, C = 0.0;
       double D = 0.0, p = 0.0, q = 0.0;
       bj = af::shared<double>(fo.size());
       af::shared<double> wj = af::shared<double>(fo.size());
       for(std::size_t i=0; i<fo.size(); i++) {
         wj[i] = 1.0;
         if(int(cf[i]) == 0) wj[i] = 2.0;
         SUMwj += wj[i];
       }
       for(std::size_t i=0; i<fo.size(); i++) {
         double fo_fo_eps = fo[i] * fo[i] / eps[i];
         double fm_fm_eps = fm[i] * fm[i] / eps[i];
         double fo_fm_eps = fm[i] * fo[i] / eps[i];
         A += ( wj[i] * fm_fm_eps );
         B += ( wj[i] * fo_fo_eps );
         C += ( wj[i] * fo_fm_eps );
         D += ( wj[i] * fo_fm_eps * fo_fm_eps );
         p += ( wj[i] * fm_fm_eps * fm_fm_eps );
         q += ( wj[i] * fo_fo_eps * fo_fo_eps );
         bj[i] = fo_fm_eps;
       }
       A /= SUMwj;
       B /= SUMwj;
       C /= SUMwj;
       D /= SUMwj;
       p /= SUMwj;
       q /= SUMwj;
       double r = (p - A * A) * (q - B * B);
       //MMTBX_ASSERT( r > 0.0 );
       if(r <= 0.0) OMEGAi = 0.0;
       else OMEGAi = (D - A * B) / std::sqrt(r);
       wi = A * B - C * C;
     }

  af::shared<double> smooth(af::shared<double> x)
     {
       if (x.size() > 1) {
         double topt1 = x[0];
         double topt2 = x[1];
         for(std::size_t i=1; i<x.size()-1; i++) {
            double topt3 = x[i+1];
            x[i]  = (topt1+topt2+topt3)/3.0;
            topt1 = topt2;
            topt2 = topt3;
         }
       }
       return x;
     }

  void alpha_beta_in_zones(af::shared<double> A_in_zones,
                           af::shared<double> B_in_zones,
                           af::shared<double> topt)
     {
       for(std::size_t i=0; i<A_in_zones.size(); i++) {
          double Ai = A_in_zones[i];
          double Bi = B_in_zones[i];
          if (topt[i] == 0.0) {
             alpha_[i] = 0.0;
             beta_[i]  = Bi;
          }
          else if (topt[i] >= 1.0e+10) {
             alpha_[i] = std::sqrt(Ai/Bi);
             beta_[i]  = 1.E-10;
          }
          else {
             double tt = 2.0 * topt[i];
             double ww = std::sqrt(1.0 + Ai * Bi * tt * tt);
             double hbeta = Bi / (ww + 1.0);
             alpha_[i] = std::sqrt(hbeta * (ww-1.0) / Ai);
             beta_[i]  = 2.0 * hbeta;
          }
       }
     }

  double solvm()
     {
       // if remove "*" and "/" below, alpha & beta will be wrong for 3 models
       // per about 1000 models
       int nst1 = 10*5;
       int nst2 = 20*5;
       double eps = 0.0001/10.;
       // default if no solution found by the iterative procedure (top, fgtop)
       double top = 0., fgtop = 0.;
       // tl is a point where the function is positive
       double tl=C/wi;
       double fgtl = funcgm(tl);
       double fgtr, tr;
       // now  funcgm(tl)>0
       for(std::size_t n1 = 0; n1 < nst1; n1++) {
           // search for the point ("left") where funcf<0
           tr=tl;
           fgtr=fgtl;
           tl=tr*0.5;
           fgtl=funcgm(tl);
           if (fgtl == 0.) {
              // exact solution found
               top=tl;
               fgtop=0.;
               return top;
           }
           else if (fgtl < 0.) {
                    // found conditions  funcg(tl)<0  and  funcg(tr)>0
                    // search for an approximate solution inside this interval
                    top=tr;
                    fgtop=fgtr;
                    for(std::size_t n2 = 0; n2 < nst2; n2++) {
                        if ((tr-tl) < eps*top) return top;
                        top=(tl*fgtr-tr*fgtl)/(fgtr-fgtl);
                        fgtop=funcgm(top);
                        if (fgtop > 0.) {
                          tr=top;
                          fgtr=fgtop;
                        }
                        else {
                          tl=top;
                          fgtl=fgtop;
                        }
                    }
           }
       }
       return top;
     }

  double funcgm(double t)
     {
       return std::sqrt(1.+4.*A*B*t*t)-1.-2.*t*blamm(t);
     }

  double blamm(double t)
     {
       double s = 0.0;
       for (std::size_t i = 0; i<bj.size(); i++) {
         if (cf[i] != 0) s += 1.0*bj[i]*fom(t*bj[i], cf[i]);
         else            s += 2.0*bj[i]*fom(t*bj[i], cf[i]);
       }
       return s/SUMwj;
     }

   double fom(double t, int lcent)
     {
       if (lcent == 0) return scitbx::math::bessel::i1_over_i0(2.*t);
       else            return std::tanh(t);
     }


  af::shared<double> alpha() { return alpha_; };
  af::shared<double> beta() { return beta_; };
protected:
  af::shared<double> fo, fm, epsilons;
  af::shared<cctbx::miller::index<> > hkl;
  double SUMwj,A,B,C,OMEGAi,wi;
  af::shared<double> alpha_,beta_,bj;
  af::shared<double> eps;
  af::shared<bool> cf;
};


class f_star_w_star_mu_nu_one_h {
public:
  f_star_w_star_mu_nu_one_h(double fo,
                            double fm,
                            double alpha,
                            double beta,
                            int    eps,
                            bool   cf_)
     {
       cf  = cf_;
       MMTBX_ASSERT(fo > 0. && fm > 0. && alpha > 0. && beta > 0.);
       MMTBX_ASSERT(eps > 0. && (cf == 0 || cf == 1));
       double w = std::sqrt(eps * beta);
       double p = fo / w;
       mu_nu(p);
       f_star_ = w * mu_ / alpha;
       w_star_ = (alpha / w) * (alpha / w) * nu_;
       if (cf == 0) w_star_ *= 2.0;
       if (f_star_ < 1.0e-6) f_star_ = 0.0;
     }

    void mu_nu(double p_)
      {
    //....................................................................
    // V.Lunin, P.Afonine & A.Urzhumtsev. Acta Cryst. (2002). A58, 270-282
    //....................................................................
       double p = std::abs(p_);
       double xstr = 0.0;
       //  ACENTRIC
       if (cf == 0) {
           if (p <= 1.0) {
              mu_ = 0.0;
              nu_ = 1.0-p*p;
           }
           else {
              if (p <= 1.3) xstr = app_as5(p);
              else          xstr = app_al4(p);
              mu_ = ref_mu(p, xstr);
              nu_ = 2.0 * (1.0 - p * p + mu_ * mu_);
           }
       }
       //  CENTRIC
       if (cf != 0) {
          if (p <= 1.0) mu_ = 0.0;
          else {
             if (p <= 1.3) xstr = app_cs5(p);
             else          xstr = app_cl4(p);
             mu_ = ref_mu(p, xstr);
          }
          nu_ = 1.0 - p * p + mu_ * mu_;
       }
     }

  double ref_mu(double p, double xstr)
     {
       double arg = p * xstr;
       double tnh = fom(arg, cf);
       if (cf == 0)
          return xstr-(p*tnh-xstr)/(2.0*p*p*(1.0-tnh/(2.0*arg)-tnh*tnh)-1.0);
       else
          return xstr-(p*tnh-xstr)/(p*p*(1.0-tnh*tnh)-1.0);
     }

  double app_cs5(double p)
     {
       double A[] = {-0.55, 0.6944643, -0.6409018, 0.6372297};
       MMTBX_ASSERT( p >= 1 );
       double p1 = p-1.0;
       double b = std::sqrt(6.0*p1);
       return b*(1+A[0]*p1+A[1]*p1*p1+A[2]*p1*p1*p1+A[3]*p1*p1*p1*p1);
     }

  double app_cl4(double p) {
       double A[] = {-2.0,2.0,-8.0,-2.0,24.0,-48.0,2.0,-48.0,256.0,-341.3333};
       MMTBX_ASSERT( p >= 1 );
       double s  = p*p;
       double es = std::exp(-2.0*s);
       return p*(1.0+A[0]*es+(A[1]+A[2]*s)*es*es+(A[3]+A[4]*s+A[5]*s*s)*es*es*es
              +(A[6]+A[7]*s+A[8]*s*s+A[9]*s*s*s)*es*es*es*es);
     }

  double app_as5(double p) {
       double A[] = {2.0,-0.8333333, 1.381944, -1.231597, 1.126676};
       MMTBX_ASSERT( p >= 1 );
       double p1=p-1.0;
       double b = std::sqrt(p1);
       return b*(A[0]+A[1]*p1+A[2]*p1*p1+A[3]*p1*p1*p1+A[4]*p1*p1*p1*p1);
     }

  double app_al4(double p) {
    double A[] = {-0.25,-0.09375,-0.0703125,-0.06884766};
    MMTBX_ASSERT( p >= 1 );
    double s=1.0/(p*p);
    return p*(1.+A[0]*s+A[1]*s*s+A[2]*s*s*s+A[3]*s*s*s*s);
     }

  double fom(double t, int lcent)
   // calculates the ratio I1(2t)/I0(2t) if lcent=0 and th(t) in other cases
   //  I1() - Modified Bessel Function of 1 order
   //  I0() - Modified Bessel Function of 0 order
   //
     {
       if (lcent == 0) return scitbx::math::bessel::i1_over_i0(2.*t);
       else            return std::tanh(t);
     }

  double w_star_one_h() { return w_star_; };
  double f_star_one_h() { return f_star_; };
  double mu_one_h()     { return mu_; };
  double nu_one_h()     { return nu_; };
protected:
  double w_star_, f_star_, mu_, nu_;
  bool  cf;
};


class f_star_w_star_mu_nu {
public:
  f_star_w_star_mu_nu(af::const_ref<double> const& fo,
                      af::const_ref<double> const& fm,
                      af::const_ref<double> const& alpha,
                      af::const_ref<double> const& beta,
                      cctbx::sgtbx::space_group const& sg,
                      af::const_ref<cctbx::miller::index<> > const& hkl)
     {
       MMTBX_ASSERT(fo.size() > 0);
       MMTBX_ASSERT(fo.size() == fm.size());
       MMTBX_ASSERT(alpha.size() == beta.size());
       MMTBX_ASSERT(fo.size() == alpha.size());
       MMTBX_ASSERT(fo.size() == hkl.size());
       eps = sg.epsilon(hkl);
       cf = sg.is_centric(hkl);
       f_star_ = af::shared<double> (fo.size());
       w_star_ = af::shared<double> (fo.size());
       mu_     = af::shared<double> (fo.size());
       nu_     = af::shared<double> (fo.size());
       nzero = 0;
       for(std::size_t i=0; i < fo.size(); i++) {
         f_star_w_star_mu_nu_one_h obj(fo[i],fm[i],alpha[i],beta[i],eps[i],cf[i]);
         f_star_[i] = obj.f_star_one_h();
         w_star_[i] = obj.w_star_one_h();
         mu_[i]     = obj.mu_one_h();
         nu_[i]     = obj.nu_one_h();
         if(f_star_[i] < 1.e-6) nzero +=1;
       }
       double w_star_max = af::max(w_star_);
       for(std::size_t i=0; i < fo.size(); i++) {
         w_star_[i] = w_star_[i] / w_star_max;
       }
     }

  af::shared<double> w_star() { return w_star_; };
  af::shared<double> f_star() { return f_star_; };
  int number_of_f_star_zero() { return nzero; };
  double ls_star_1()          { return ls_star_1_; };
  double ls_star_2()          { return ls_star_2_; };
  double k_1()                { return k_1_; };
  double k_2()                { return k_2_; };
  af::shared<double> mu()     { return mu_; };
  af::shared<double> nu()     { return nu_; };
protected:
  af::shared<double> w_star_, f_star_, fom_, mu_, nu_;
  af::shared<int> eps;
  af::shared<bool> cf;
  double ls_star_1_, ls_star_2_, k_1_, k_2_;
  int nzero;
};

class fom_and_phase_error {
public:
  fom_and_phase_error(af::shared<double> const& fo,
                      af::shared<double> const& fm,
                      af::shared<double> const& alpha,
                      af::shared<double> const& beta,
                      af::shared<double> const& epsilons,
                      af::shared<bool> const& centric_flags
                      )
     {
       fo_ = fo;
       fm_ = fm;
       alpha_ = alpha;
       beta_ = beta;
       MMTBX_ASSERT( fo.size() > 0 && (fo.size() == fm.size()) );
       MMTBX_ASSERT( alpha.size() == beta.size() );
       for(std::size_t i=0; i < epsilons.size(); i++) {
         eps_.push_back(epsilons[i]);
         cf_.push_back(centric_flags[i]);
       }
     }

  //! calculate figures of merit (fom)
  //! formula 6 and 9 in Lunin & Skovoroda (Acta Cryst. (1995). A51, 880-887)
  af::shared<double> fom()
     {
       fom_ = af::shared<double> (fo_.size());
       for(std::size_t i=0; i < fo_.size(); i++) {
         if(beta_[i] > 0.0 && alpha_[i] >= 0.0) {
            double p = alpha_[i] * fo_[i] * fm_[i] / (eps_[i] * beta_[i]);
            if (cf_[i] == 0) fom_[i] = scitbx::math::bessel::i1_over_i0(2.*p);
            else             fom_[i] = std::tanh(p);
         }
         else {
            fom_[i] = 0.0;
         }
       }
       return fom_;
     }

  //! calculate absolute phase error <|phi - phi_model|>
  //! formula 7 and 10 in Lunin & Skovoroda (Acta Cryst. (1995). A51, 880-887)
  af::shared<double> phase_error()
     {
       double mp = 180./(scitbx::constants::pi*scitbx::constants::pi);
       phase_error_ = af::shared<double> (fo_.size());
       for(std::size_t i=0; i < fo_.size(); i++) {
         if(beta_[i] > 0.0 && alpha_[i] >= 0.0) {
            double p = alpha_[i] * fo_[i] * fm_[i] / (eps_[i] * beta_[i]);
            if (cf_[i] == 0) {
                phase_error_[i] = mp * simp(2.*p) / pseudo_i0(2.*p);
            }
            else {
                phase_error_[i] = 180.0 * std::exp(-2.*p) / (1.0 + std::exp(-2.*p));
            }
         }
         else {
            phase_error_[i] = 90.0;
         }
       }
       return phase_error_;
     }

  static
  double fint(double phi, double tau)
     {
      if (tau <= 3.75) return phi*std::exp(tau*std::cos(phi));
      else             return phi*std::exp(tau*(std::cos(phi)-1.));
     }

  //! calculate I0(x): very specific to this task. Avoids overflow problem.
  //! return I0(x) if |x| < 3.75
  //! return exp(-x)*I0(x) if x >= 3.75 (*** not I0(x) ***)
  //! "missing" exponent appears in fint() function
  static
  double pseudo_i0(double const& x)
     {
       double abs_x = x;
       if (abs_x < 0) abs_x = -abs_x;
       double bessel_i0;
       if (abs_x/3.75 < 1.0) {
          double y = x/3.75;
          y *= y;
          bessel_i0=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+
                    y*(0.2659732+y*(0.0360768+y*0.0045813)))));
       }
       else {
          double y = 3.75/abs_x;
          y=0.39894228+y*(0.01328592+y*(0.00225319+y*(-0.00157565+
            y*(0.00916281+y*(-0.02057706+y*(0.02635537+
            y*(-0.01647633+y*0.00392377)))))));
          bessel_i0=y/std::sqrt(abs_x);
       }
       return bessel_i0;
     }

  //! integral by Simpson method
  static
  double simp(double tau)
     {
       double hhh,hhhh,x4,x5,f4,f5,q2,q3,qqq;
       af::shared<double> q(30),f2s(30),f5s(30),x2s(30),x5s(30);
       double a = 0.0, b = 3.1415926, t = 0.0001, t1 = 0.00001;
       double tt=t*t;
       double h=(b-a)/2.;
       int i=0;
       double qq=0.;
       double x1=a;
       double x2=b;
       double f1 = fint(x1,tau);
       double f2 = fint(x2,tau);
       double x3 = x1+h;
       double f3 = fint(x3,tau);
       double q1 = (f1+f2+4.*f3)*h;
       double eps=t*std::abs(q1);
       if(std::abs(q1)<=t) eps=tt;
       lbl_2:;
       h=h/2.;
       hhh=std::abs(x1);
       hhhh=std::abs(x2);
       if(hhhh>hhh) hhh=hhhh;
       if(std::abs(h)/hhh<=t1) goto lbl_4;
       x4=x1+h;
       x5=x3+h;
       f4 = fint(x4,tau);
       f5 = fint(x5,tau);
       q2=(f1+f3+4.*f4)*h;
       q3=(f3+f2+4.*f5)*h;
       qqq=q2+q3;
       if(std::abs(qqq-q1)<=eps) goto lbl_1;
       if(i>=30) goto lbl_1;
       i=i+1;
       f2s[i]=f2;
       x2s[i]=x2;
       f5s[i]=f5;
       x5s[i]=x5;
       q[i]=q3;
       f2=f3;
       x2=x3;
       f3=f4;
       x3=x4;
       q1=q2;
       goto lbl_2;
       lbl_1:;
       qq=qq+qqq;
       eps=t*std::abs(qq);
       if(std::abs(qq)<=t) eps=tt;
       lbl_4:;
       if(i==0) goto lbl_3;
       x1=x2;
       f1=f2;
       x2=x2s[i];
       f2=f2s[i];
       x3=x5s[i];
       f3=f5s[i];
       q1=q[i];
       i=i-1;
       h=(x2-x1)/2.;
       goto lbl_2;
       lbl_3:;
       return qq/3.;
     }

protected:
  af::shared<double> fo_, fm_, alpha_, beta_;
  af::shared<double> fom_, phase_error_;
  af::shared<double> eps_;
  af::shared<bool> cf_;
};

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
class sasha_error_calculator{
public:
      sasha_error_calculator(af::const_ref<vec3<double> > const& r1f,
                             af::const_ref<vec3<double> > const& r1c,
                             af::const_ref<vec3<double> > const& r2f,
                             af::const_ref<vec3<double> > const& r2c,
                             af::const_ref<std::string> const& lab1,
                             af::const_ref<std::string> const& lab2,
                             cctbx::uctbx::unit_cell const& uc,
                             cctbx::sgtbx::space_group const& sg,
                             double rad=3.0)
      {
       int nat1 = r1f.size();
       int nat2 = r2f.size();
       MMTBX_ASSERT(nat1 >= nat2);
       MMTBX_ASSERT(r1f.size()==r1c.size() && r1c.size()==lab1.size());
       MMTBX_ASSERT(r2f.size()==r2c.size() && r2c.size()==lab2.size());
       //
       //     one-to-one identification of atoms in two files
       //     nrefa1: atom from the 1st file to this atom in the 2nd file
       //     nrefa2: atom from the 2nd file to this atom in the 1st file
       //     9999 stands for the removed atom
       //
       af::shared<int> nrefa1 = af::shared<int> (nat1);
       af::shared<int> nrefa2 = af::shared<int> (nat1);
       int j1 = 0;
       for(int i2=0; i2 < nat2; i2++) {
           std::string title = lab2[i2];
           for(int i1=j1; i1 < nat1; i1++) {
               if(lab1[i1] == title) {
                  nrefa2[i2]=i1;
                  nrefa1[i1]=i2;
                  j1=i1+1;
                  break;
               }
               else {
                  nrefa1[i1] = 9999999;
               }
           }
       }
       if(j1 <= nat1) {
          for(std::size_t i1=j1; i1 < nat1; i1++) {
              nrefa1[i1] = 9999999;
          }
       }
       //
       // expansion of the 1st model by symmetries and translations
       // to generate all atoms in the unit cell and nearby
       //
       // ATTENSION : this is a time-consuming part.
       // if several models with be checked during the same refinement run
       // probably it is better to save in a file the information
       // from this block and read it for other runs...
       //
       // region around the unit cell where interacting atoms may be
       // PA  orthogonal cell is considered
       //
       double xmin=  -rad/uc.parameters()[0];
       double xmax=1.+rad/uc.parameters()[0];
       double ymin=  -rad/uc.parameters()[1];
       double ymax=1.+rad/uc.parameters()[1];
       double zmin=  -rad/uc.parameters()[3];
       double zmax=1.+rad/uc.parameters()[3];
       af::shared<double> xc1;
       af::shared<double> yc1;
       af::shared<double> zc1;
       af::shared<int> numc1;
       // symmetries:
       double rs[12][240];
       using cctbx::sgtbx::rt_mx;
       using cctbx::sgtbx::rot_mx;
       using cctbx::sgtbx::tr_vec;
       for (int j = 0; j < sg.order_z(); j++) {
            rt_mx rt = sg(j);
            rot_mx r = rt.r();
            tr_vec t = rt.t();
            for(int i=0; i<9; i++) {
                rs[i][j] = r.num()[i] / static_cast<double>(r.den());
            }
            for(int i=0; i<3; i++) {
                rs[i+9][j] = t.num()[i] / static_cast<double>(t.den());
            }
       }
       // checking the atoms
       int jc=0;
       for(int i1=0; i1 < nat1; i1++) {
           double xf = r1f[i1][0];
           double yf = r1f[i1][1];
           double zf = r1f[i1][2];
           // reduction of cryst. (fractional) coords. to the unit cell
           // now x,y,z vary between 0 and 1
           int ix=int(xf);
           if(xf < 0.0) ix=ix-1;
           xf=xf-ix;
           int iy=int(yf);
           if(yf < 0.0) iy=iy-1;
           yf=yf-iy;
           int iz=int(zf);
           if(zf < 0.0) iz=iz-1;
           zf=zf-iz;
           // atom expansion over symmetries
           for (int k = 0; k < sg.order_z(); k++) {
               double xs = rs[0][k]*xf+rs[1][k]*yf+rs[2][k]*zf+rs[9][k];
               double ys = rs[3][k]*xf+rs[4][k]*yf+rs[5][k]*zf+rs[10][k];
               double zs = rs[6][k]*xf+rs[7][k]*yf+rs[8][k]*zf+rs[11][k];
               // new point may be out of the cell - correct if necessary
               if(xs <  0.) xs=xs+1.;
               if(xs >= 1.) xs=xs-1.;
               if(ys <  0.) ys=ys+1.;
               if(ys >= 1.) ys=ys-1.;
               if(zs <  0.) zs=zs+1.;
               if(zs >= 1.) zs=zs-1.;
               // shift over neighbour cells
               for (int iz = -1; iz <=1; iz++) {
                 double z0=zs+iz;
                 if(z0>zmin && z0<zmax) {
                   for (int iy = -1; iy <=1; iy++) {
                     double y0=ys+iy;
                     if(y0>ymin && y0<ymax) {
                       for (int ix = -1; ix <=1; ix++) {
                         double x0=xs+ix;
                         if(x0>xmin && x0<xmax) {
                           jc=jc+1;
                           xc1.push_back(x0);
                           yc1.push_back(y0);
                           zc1.push_back(z0);
                           numc1.push_back(i1);
                         }
                       }
                     }
                   }
                 }
               }
           }
       }
       int natex1=jc;
       // orthogonalization of the obtained coordinates
       for (int i = 0; i <natex1; i++) {
            vec3<double> const& vo = uc.orthogonalize(
              vec3<double>(xc1[i],yc1[i],zc1[i]));
            xc1[i] = vo[0];
            yc1[i] = vo[1];
            zc1[i] = vo[2];
       }
       // reduction of atoms of the second model into the basic unit cell
       //
       // checking the atoms
       af::shared<double> xc2;
       af::shared<double> yc2;
       af::shared<double> zc2;
       for (int i2 = 0; i2 <nat2; i2++) {
           double xf = r2f[i2][0];
           double yf = r2f[i2][1];
           double zf = r2f[i2][2];
           //reduction of cryst. (fractional) coords. to the unit cell
           //now x,y,z vary between 0 and 1
           int ix=int(xf);
           if(xf<0.0) ix=ix-1;
           xf=xf-ix;
           int iy=int(yf);
           if(yf<0.0) iy=iy-1;
           yf=yf-iy;
           int iz=int(zf);
           if(zf<0.0) iz=iz-1;
           zf=zf-iz;
           // coordinate orthogonalisation
           vec3<double> const& vo = uc.orthogonalize(vec3<double>(xf,yf,zf));
           xc2.push_back(vo[0]);
           yc2.push_back(vo[1]);
           zc2.push_back(vo[2]);
       }
       // computing the usual errors
       af::shared<double> adist;
       for (int i2 = 0; i2 <nat2; i2++) {
         int i1 = nrefa2[i2];
         double dx=r1c[i1][0]-r2c[i2][0];
         double dy=r1c[i1][1]-r2c[i2][1];
         double dz=r1c[i1][2]-r2c[i2][2];
         double dist = std::sqrt(dx*dx+dy*dy+dz*dz);
         adist.push_back(dist);
       }
       // computing the optimal errors
       double rad2=rad*rad;
       af::shared<double> adopt;
       for (int i2 = 0; i2 <nat2; i2++) {
          double xx=xc2[i2];
          double yy=yc2[i2];
          double zz=zc2[i2];
          double dd=rad2;
       // checking the atoms from the 1st expandend model
          for (int i1 = 0; i1 <natex1; i1++) {
              double dx=xc1[i1]-xx;
              double dy=yc1[i1]-yy;
              double dz=zc1[i1]-zz;
              if(dx>-rad && dx<rad && dy>-rad && dy<rad && dz>-rad && dz<rad) {
                 double dist2 = dx*dx+dy*dy+dz*dz;
                 if(dist2<dd) {
                   dd=dist2;
                 }
              }
          }
          adopt.push_back(std::sqrt(dd));
       }
       // get statistics of the error distribution
       distm_=0;
       doptm_=0;
       for (int i2 = 0; i2 <nat2; i2++) {
         distm_=distm_+adist[i2];
         doptm_=doptm_+adopt[i2];
       }
       distm_=distm_/nat2;
       doptm_=doptm_/nat2;

      } //end constructor

      double doptm() { return doptm_; };
      double distm() { return distm_; };

protected:
      double distm_,doptm_;
};
///////////////////////////////////////////////////////////
class peak_clustering{
public:
    peak_clustering(af::const_ref<vec3<double> > const& r1f,
                    af::const_ref<vec3<double> > const& r2f,
                    af::const_ref<double> const& h1,
                    af::const_ref<double> const& h2,
                    cctbx::uctbx::unit_cell const& uc,
                    double const& cutoff)
    {
      for(int i = 0; i < r1f.size(); i++) {
          cctbx::fractional<double> s = vec3<double>(0,0,0);
          double h = 0.0;
          for(int j = 0; j < r2f.size(); j++) {
              cctbx::fractional<double> site_1 = r1f[i];
              cctbx::fractional<double> site_2 = r2f[j];
              double d = uc.distance(site_1, site_2);
              if(d <= cutoff) {
                 s += (site_1 * h1[i] + site_2 * h2[j]);
                 h += (h1[i] + h2[j]);
              }
          }
          if(h > 0.0) {
             sites_.push_back(s/h);
             heights_.push_back(h/s.size());
          }
      }
    }
    af::shared<vec3<double> > sites() { return sites_; };
    af::shared<double> heights() { return heights_; };

protected:
    af::shared<vec3<double> > sites_;
    af::shared<double> heights_;
};

///////////////////////////////////////////////////////////

af::shared<double> fo_fc_alpha_over_eps_beta(
                     af::shared<double> const& fo,
                     af::shared<double> const& fm,
                     af::shared<double> const& alpha,
                     af::shared<double> const& beta,
                     cctbx::sgtbx::space_group const& sg,
                     af::const_ref<cctbx::miller::index<> > hkl);

}} // namespace max_lik
#endif
