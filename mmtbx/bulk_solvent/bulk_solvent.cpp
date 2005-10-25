#include <assert.h>
#include <math.h>
#include <iostream>
#include <mmtbx/bulk_solvent/bulk_solvent.h>
#include <mmtbx/error.h>
#include <cmath>
#include <scitbx/math/bessel.h>
#include <cctbx/xray/targets.h>
#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>


using namespace std;
namespace mmtbx { namespace bulk_solvent {

vec3<double> ksol_bsol_grid_search(
                             af::const_ref<double> const& fo,
                             af::const_ref< std::complex<double> > const& fc,
                             af::const_ref< std::complex<double> > const& fm,
                             sym_mat3<double> const& u,
                             af::const_ref<double> const& ksol_range,
                             af::const_ref<double> const& bsol_range,
                             double const& r_ref,
                             af::const_ref<cctbx::miller::index<> > const& hkl,
                             cctbx::uctbx::unit_cell const& uc)
{
 MMTBX_ASSERT(hkl.size() == fo.size());
 MMTBX_ASSERT(fo.size() == fc.size() && fc.size() == fm.size());
 double k_best = 0.0;
 double b_best = 0.0;
 double r_best = r_ref;
 for(std::size_t i=0; i < ksol_range.size(); i++) {
   for(std::size_t j=0; j < bsol_range.size(); j++) {
     double r = r_factor_aniso_fast(fo,fc,fm,u,ksol_range[i],bsol_range[j],hkl,uc);
     if(r < r_best) {
       k_best = ksol_range[i];
       b_best = bsol_range[j];
       r_best = r;
     }
   }
 }
 return vec3<double> (k_best,b_best,r_best);
}

double fu_star(sym_mat3<double> const& u,
               cctbx::miller::index<> const& mi)
{
    double arg = -0.25 * (u[0]*mi[0]*mi[0] +
                          u[1]*mi[1]*mi[1] +
                          u[2]*mi[2]*mi[2] +
                          u[3]*mi[0]*mi[1] +
                          u[4]*mi[0]*mi[2] +
                          u[5]*mi[1]*mi[2]);
    if(arg > std::log(5.0)) return 1.0;
    return std::exp(arg);
}

af::shared<double> fu_aniso(sym_mat3<double> const& u,
                            af::const_ref<cctbx::miller::index<> > const& hkl,
                            cctbx::uctbx::unit_cell const& uc)
{
    mat3<double> a = uc.fractionalization_matrix();
    sym_mat3<double> u_star = sym_mat3<double> (u).tensor_transform(a);
    af::shared<double> fu_mem(hkl.size(), af::init_functor_null<double>());
    double* fu = fu_mem.begin();
    for(std::size_t i=0; i < hkl.size(); i++) {
      fu[i] = fu_star(u_star, hkl[i]);
    }
    return fu_mem;
}

target_gradients_aniso::target_gradients_aniso(
                               af::const_ref<double> const& fo,
                               af::const_ref< std::complex<double> > const& fc,
                               af::const_ref< std::complex<double> > const& fm,
                               sym_mat3<double> const& u,
                               double const& ksol,
                               double const& bsol,
                               af::const_ref<cctbx::miller::index<> > const& hkl,
                               cctbx::uctbx::unit_cell const& uc,
                               bool const& calc_grad_u,
                               bool const& calc_grad_ksol,
                               bool const& calc_grad_bsol,
                               bool const& trace_zero)
{
    af::shared<double> s_mem = uc.stol_sq(hkl);
    af::const_ref<double> s = s_mem.const_ref();
    mat3<double> a = uc.fractionalization_matrix();
    sym_mat3<double> u_star = u.tensor_transform(a);
    MMTBX_ASSERT(hkl.size()==fo.size());
    MMTBX_ASSERT(fo.size()==fc.size()&&fc.size()==fm.size());
    std::vector<double> f_b(fo.size());
    std::vector<double> fmodel_complex_abs(fo.size());
    double num=0.0;
    double denum=0.0;
    for(std::size_t i=0; i < fo.size(); i++) {
      double fu = fu_star(u_star, hkl[i]);
      f_b[i]=std::exp(-bsol * s[i]);
      double fmodel_abs = std::abs((fc[i] + f_b[i]*ksol * fm[i]) * fu);
      fmodel_complex_abs[i]=fmodel_abs;
      num += fo[i] * fmodel_abs;
      denum += fmodel_abs*fmodel_abs;
    }
    MMTBX_ASSERT(denum > 0.0);
    double sc = num/denum;
    tgx = 0.0;
    gtgx_ksol = 0.0;
    gtgx_bsol = 0.0;
    scale_tgx = 0.0;
    MMTBX_ASSERT(gtgx_u.size() == 0);
    gtgx_u.resize(6, 0.);
    double* gxu = gtgx_u.begin();
    for(std::size_t i=0; i < fo.size(); i++) {
      double fmodel = fmodel_complex_abs[i];
      scale_tgx += fo[i]*fo[i];
      double delta = fo[i] - fmodel * sc;
      tgx += delta*delta;
//grads:
      if(calc_grad_u) {
        double coeff = -delta * 2 * fmodel * sc * (-0.25);
        cctbx::miller::index<> const& mi = hkl[i];
        double t1 = a[0]*mi[0] + a[3]*mi[1] + a[6]*mi[2];
        double t2 = a[1]*mi[0] + a[4]*mi[1] + a[7]*mi[2];
        double t3 = a[2]*mi[0] + a[5]*mi[1] + a[8]*mi[2];
        gxu[0] += coeff * t1*t1;
        gxu[1] += coeff * t2*t2;
        gxu[2] += coeff * t3*t3;
        gxu[3] += coeff * 2.*t1*t2;
        gxu[4] += coeff * 2.*t1*t3;
        gxu[5] += coeff * 2.*t2*t3;
      }
      if(calc_grad_ksol || calc_grad_bsol) {
        double uvs_plus_usv = std::real(fc[i]*std::conj(fm[i])+fm[i]*std::conj(fc[i]));
        double mod_v_sq = std::abs(fm[i]) * std::abs(fm[i]);
        double coeff = (uvs_plus_usv + f_b[i]*ksol * mod_v_sq * 2) / (fmodel * 2);
        if(calc_grad_ksol) { gtgx_ksol +=-delta * 2 * coeff * f_b[i]; }
        if(calc_grad_bsol) { gtgx_bsol += delta * 2 * coeff * f_b[i]*ksol * s[i]; }
      }
   }
   tgx /= scale_tgx;
   gtgx_ksol /= scale_tgx;
   gtgx_bsol /= scale_tgx;
   for(std::size_t j=0; j < 6; j++) {
     gxu[j] /= scale_tgx;
   }
   if(trace_zero) {
      double trace_zero_tg  = (u[0]+u[1]+u[2]) * (u[0]+u[1]+u[2]);
      double trace_zero_gtg = (u[0]+u[1]+u[2]) * 2.0;
      tgx += trace_zero_tg;
      gxu[0] += trace_zero_gtg;
      gxu[1] += trace_zero_gtg;
      gxu[2] += trace_zero_gtg;
   }
}

double scale_factor_aniso(af::const_ref<double> const& fo,
                          af::const_ref< std::complex<double> > const& fc,
                          af::const_ref< std::complex<double> > const& fm,
                          sym_mat3<double> const& u,
                          double const& ksol,
                          double const& bsol,
                          af::const_ref<cctbx::miller::index<> > const& hkl,
                          cctbx::uctbx::unit_cell const& uc,
                          af::const_ref<double> const& s)
{
    MMTBX_ASSERT(fo.size()==fc.size()&&fc.size()==fm.size());
    MMTBX_ASSERT(hkl.size()==fo.size());
    mat3<double> a = uc.fractionalization_matrix();
    sym_mat3<double> u_star = sym_mat3<double> (u).tensor_transform(a);
    double num =0.0;
    double denum=0.0;
    for(std::size_t i=0; i < fo.size(); i++) {
      double f_u = fu_star(u_star, hkl[i]);
      double func_k_b = ksol * std::exp(-bsol * s[i]);
      double fmodel = std::abs(fc[i] + func_k_b * fm[i]) * f_u;
      num += fo[i] * fmodel;
      denum += fmodel*fmodel;
    }
    MMTBX_ASSERT(denum > 0.0);
    return num/denum;
}

double r_factor_aniso(af::const_ref<double> const& fo,
                      af::const_ref< std::complex<double> > const& fc,
                      af::const_ref< std::complex<double> > const& fm,
                      sym_mat3<double> const& u,
                      double const& ksol,
                      double const& bsol,
                      af::const_ref<cctbx::miller::index<> > const& hkl,
                      cctbx::uctbx::unit_cell const& uc,
                      double const& sc,
                      af::const_ref<double> const& s)
{
    MMTBX_ASSERT(hkl.size()==fo.size());
    MMTBX_ASSERT(fo.size()==fc.size()&&fc.size()==fm.size());
    MMTBX_ASSERT(u.size()==6);
    mat3<double> a = uc.fractionalization_matrix();
    sym_mat3<double> u_star = sym_mat3<double> (u).tensor_transform(a);
    double num =0.0;
    double denum=0.0;
    for(std::size_t i=0; i < fo.size(); i++) {
      double f_u = fu_star(u_star, hkl[i]);
      double func_k_b = ksol * std::exp(-bsol * s[i]);
      double fmodel = std::abs(fc[i] + func_k_b * fm[i]) * f_u;
      num += std::abs(fo[i] - fmodel * sc);
      denum += fo[i];
    }
    MMTBX_ASSERT(denum > 0.0);
    return num/denum;
}

double r_factor_aniso_fast(af::const_ref<double> const& fo,
                           af::const_ref< std::complex<double> > const& fc,
                           af::const_ref< std::complex<double> > const& fm,
                           sym_mat3<double> const& u,
                           double const& ksol,
                           double const& bsol,
                           af::const_ref<cctbx::miller::index<> > const& hkl,
                           cctbx::uctbx::unit_cell const& uc)
{
    af::shared<double> s_mem = uc.stol_sq(hkl);
    af::const_ref<double> s = s_mem.const_ref();
    MMTBX_ASSERT(hkl.size()==fo.size());
    MMTBX_ASSERT(fo.size()==fc.size()&&fc.size()==fm.size());
    mat3<double> a = uc.fractionalization_matrix();
    sym_mat3<double> u_star = sym_mat3<double> (u).tensor_transform(a);
    std::vector<double> fmodel_abs(fo.size());
    double num=0.0;
    double denum=0.0;
    for(std::size_t i=0; i < fo.size(); i++) {
      double fu = fu_star(u_star, hkl[i]);
      double f_b_k = std::exp(-bsol * s[i]) * ksol;
      double fmodel_abs_ = std::abs((fc[i] + f_b_k * fm[i]) * fu);
      fmodel_abs[i]=fmodel_abs_;
      num += fo[i] * fmodel_abs_;
      denum += fmodel_abs_*fmodel_abs_;
    }
    MMTBX_ASSERT(denum > 0.0);
    double sc = num/denum;
    num =0.0;
    denum=0.0;
    for(std::size_t i=0; i < fo.size(); i++) {
      num += std::abs(fo[i] - fmodel_abs[i] * sc);
      denum += fo[i];
    }
    MMTBX_ASSERT(denum > 0.0);
    return num/denum;
}

double r_factor(af::const_ref<double> const& fo,
                af::const_ref< std::complex<double> > const& fc)
{
    MMTBX_ASSERT(fo.size()==fc.size());
    double num=0.0;
    double denum=0.0;
    std::vector<double> fc_abs(fo.size());
    for(std::size_t i=0; i < fo.size(); i++) {
      fc_abs[i]=std::abs(fc[i]);
      num += fo[i] * fc_abs[i];
      denum += fc_abs[i]*fc_abs[i];
    }
    MMTBX_ASSERT(denum > 0.0);
    double sc = num/denum;
    num =0.0;
    denum=0.0;
    for(std::size_t i=0; i < fo.size(); i++) {
      num += std::abs(fo[i] - fc_abs[i] * sc);
      denum += fo[i];
    }
    MMTBX_ASSERT(denum > 0.0);
    return num/denum;
}

double scale(af::const_ref<double> const& fo,
             af::const_ref< std::complex<double> > const& fc)
{
    MMTBX_ASSERT(fo.size()==fc.size());
    std::vector<double> fmodel_abs;
    double num=0.0;
    double denum=0.0;
    std::vector<double> fc_abs;
    for(std::size_t i=0; i < fo.size(); i++) {
      double fc_abs_ = std::abs(fc[i]);
      num += fo[i] * fc_abs_;
      denum += fc_abs_*fc_abs_;
    }
    MMTBX_ASSERT(denum > 0.0);
    return num/denum;
}

target_gradients_aniso_ml::target_gradients_aniso_ml(
                               af::const_ref<double> const& fo,
                               af::const_ref< std::complex<double> > const& fc,
                               af::const_ref< std::complex<double> > const& fm,
                               sym_mat3<double> const& u,
                               double const& ksol,
                               double const& bsol,
                               af::const_ref<cctbx::miller::index<> > const& hkl,
                               cctbx::uctbx::unit_cell const& uc,
                               cctbx::sgtbx::space_group const& sg,
                               af::const_ref<bool> const& gradient_flags,
                               af::const_ref<double> const& alpha,
                               af::const_ref<double> const& beta,
                               double k,
                               bool const& trace_zero)
{
    af::shared<double> s_mem = uc.stol_sq(hkl);
    af::const_ref<double> s = s_mem.const_ref();
    af::shared<int> eps_mem = sg.epsilon(hkl);
    af::const_ref<int> eps = eps_mem.const_ref();
    af::shared<bool> cf_mem = sg.is_centric(hkl);
    af::const_ref<bool> cf = cf_mem.const_ref();
    MMTBX_ASSERT(hkl.size()==fo.size());
    MMTBX_ASSERT(gradient_flags.size()==3 && gradient_flags[0]==gradient_flags[1]);
    bool compute_ksol_bsol_grad = gradient_flags[0];
    bool compute_uaniso_grad = gradient_flags[2];
    MMTBX_ASSERT(fo.size()==fc.size()&&fc.size()==fm.size());
    mat3<double> a = uc.fractionalization_matrix();
    sym_mat3<double> u_star = sym_mat3<double> (u).tensor_transform(a);
    std::vector<double> fu(fo.size());
    std::vector<double> f_b(fo.size());
    std::vector<double> f_b_k(fo.size());
    std::vector<double> fmodel_complex_abs(fo.size());
    for(std::size_t i=0; i < fo.size(); i++) {
      fu[i]=fu_star(u_star, hkl[i]);
      f_b[i]=std::exp(-bsol * s[i]);
      f_b_k[i]=f_b[i] * ksol;
      fmodel_complex_abs[i]=std::abs((fc[i] + f_b_k[i] * fm[i]) * fu[i]);
    }
    tgx = 0.0;
    gtgx_ksol = 0.0;
    gtgx_bsol = 0.0;
    gtgx_k = 0.0;
    MMTBX_ASSERT(gtgx_u.size() == 0);
    gtgx_u.resize(6, 0.);
    double* gxu = gtgx_u.begin();
    double coeff_1 = 0.0;
    for(std::size_t i=0; i < fo.size(); i++) {
      int cs = int(cf[i]);
      double eps_beta = eps[i]*beta[i];
      double alpha_sqr = alpha[i]*alpha[i];
      double fmodel = fmodel_complex_abs[i];
      double t_q = alpha[i]*fo[i]*fmodel/eps_beta;
      tgx += cctbx::xray::targets::maximum_likelihood_target_one_h(fo[i],fmodel,alpha[i],beta[i],k,eps[i],cs);
      if(compute_uaniso_grad || compute_ksol_bsol_grad) {
        if(fmodel > 0.0) {
           double d_p =-alpha_sqr*fmodel/eps_beta;
           double d_q = alpha[i]*fo[i]/eps_beta;
           double d_r = scitbx::math::bessel::i1_over_i0(2.*t_q);
           double d_s = std::tanh(t_q);
           coeff_1 = ( (1-cs)*(2.*d_p+2.*d_q*d_r)+cs*(d_p+d_q*d_s) )*(-1.0);
        }
      }
      if(compute_uaniso_grad) {
        double coeff_2 = coeff_1 * fmodel * (-0.25);
        cctbx::miller::index<> const& mi = hkl[i];
        double t1 = a[0]*mi[0] + a[3]*mi[1] + a[6]*mi[2];
        double t2 = a[1]*mi[0] + a[4]*mi[1] + a[7]*mi[2];
        double t3 = a[2]*mi[0] + a[5]*mi[1] + a[8]*mi[2];
        gxu[0] += coeff_2 * t1*t1;
        gxu[1] += coeff_2 * t2*t2;
        gxu[2] += coeff_2 * t3*t3;
        gxu[3] += coeff_2 * 2.*t1*t2;
        gxu[4] += coeff_2 * 2.*t1*t3;
        gxu[5] += coeff_2 * 2.*t2*t3;
      }
      if(compute_ksol_bsol_grad) {
        double uvs_plus_usv = std::real(fc[i]*std::conj(fm[i])+fm[i]*std::conj(fc[i]));
        double mod_v_sq = std::abs(fm[i]) * std::abs(fm[i]);
        double coeff_3 = (uvs_plus_usv + f_b_k[i] * mod_v_sq * 2) / (fmodel * 2);
        double coeff_1_3 = coeff_1 * coeff_3;
        gtgx_ksol +=  coeff_1_3 * f_b[i];
        gtgx_bsol += -coeff_1_3 * f_b_k[i] * s[i];
      }
      if(1) {
        gtgx_k += cctbx::xray::targets::d_maximum_likelihood_target_one_h_over_k(
                                    fo[i],fmodel,alpha[i],beta[i],k,eps[i],cs);
      }
   }
   if(trace_zero) {
      double trace_zero_tg  = (u[0]+u[1]+u[2]) * (u[0]+u[1]+u[2]);
      double trace_zero_gtg = (u[0]+u[1]+u[2]) * 2.0;
      tgx += trace_zero_tg;
      gxu[0] += trace_zero_gtg;
      gxu[1] += trace_zero_gtg;
      gxu[2] += trace_zero_gtg;
   }
}

}} // namespace mmtbx::bulk_solvent
