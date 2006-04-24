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
                             sym_mat3<double> const& b_cart,
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
     double r = r_factor_aniso_fast(fo,fc,fm,b_cart,ksol_range[i],bsol_range[j],hkl,uc);
     if(r < r_best) {
       k_best = ksol_range[i];
       b_best = bsol_range[j];
       r_best = r;
     }
   }
 }
 return vec3<double> (k_best,b_best,r_best);
}


double fu_star_(sym_mat3<double> const& b,
               cctbx::miller::index<> const& mi,
               mat3<double> a)
{
// Matamatica proven equivalent of fu_star()
    int h = mi[0];
    int k = mi[1];
    int l = mi[2];
    double qq =   a[0]*a[0]*b[0]*h*h +
  a[1]*a[1]*b[1]*h*h +
  a[2]*a[2]*b[2]*h*h +
2*a[2]*a[5]*b[2]*h*k +
2*a[2]*a[3]*b[4]*h*k +
2*a[2]*a[4]*b[5]*h*k +
  a[3]*a[3]*b[0]*k*k +
  a[4]*a[4]*b[1]*k*k +
  a[5]*a[5]*b[2]*k*k +
2*a[3]*a[4]*b[3]*k*k +
2*a[3]*a[5]*b[4]*k*k +
2*a[4]*a[5]*b[5]*k*k +
2*(a[2]*(a[8]*b[2] + a[6]*b[4] + a[7]*b[5])*h +
(a[3]*a[6]*b[0] + a[4]*a[7]*b[1] + a[5]*a[8]*b[2] +
 a[4]*a[6]*b[3] + a[3]*a[7]*b[3] + a[5]*a[6]*b[4] +
 a[3]*a[8]*b[4] + a[5]*a[7]*b[5] + a[4]*a[8]*b[5])*k)*l +
(a[7]*a[7]*b[1] + a[8]*a[8]*b[2] + a[6]*(a[6]*b[0] + 2*a[7]*b[3] + 2*a[8]*b[4]) + 2*a[7]*a[8]*b[5])*l*l +
2*a[0]*h*(a[1]*b[3]*h + a[2]*b[4]*h + a[3]*b[0]*k +
        a[4]*b[3]*k + a[5]*b[4]*k + a[6]*b[0]*l + a[7]*b[3]*l + a[8]*b[4]*l) +
2*a[1]*h*(a[2]*b[5]*h + a[4]*b[1]*k + a[3]*b[3]*k + a[5]*b[5]*k + a[7]*b[1]*l +
        a[6]*b[3]*l + a[8]*b[5]*l);
    double arg = -0.25 * (qq);
    if(arg > std::log(5.0)) return 1.0;
    return std::exp(arg);
}


double fu_star(sym_mat3<double> const& u_star,
               cctbx::miller::index<> const& mi)
{
    double arg = -0.25 * (u_star[0]*mi[0]*mi[0] +
                          u_star[1]*mi[1]*mi[1] +
                          u_star[2]*mi[2]*mi[2] +
                       2.*u_star[3]*mi[0]*mi[1] +
                       2.*u_star[4]*mi[0]*mi[2] +
                       2.*u_star[5]*mi[1]*mi[2]);
    if(arg > std::log(5.0)) return 1.0;
    return std::exp(arg);
}

af::shared<double> fb_cart(sym_mat3<double> const& b_cart,
                            af::const_ref<cctbx::miller::index<> > const& hkl,
                            cctbx::uctbx::unit_cell const& uc)
{
    mat3<double> a = uc.fractionalization_matrix();
    sym_mat3<double> u_star = sym_mat3<double> (b_cart).tensor_transform(a);
    af::shared<double> fu_mem(hkl.size(), af::init_functor_null<double>());
    double* fu = fu_mem.begin();
    for(std::size_t i=0; i < hkl.size(); i++) {
      double zz = fu_star_(b_cart, hkl[i], a);
      fu[i] = fu_star(u_star, hkl[i]);
      // Make sure both functions are equal.
      //std::cout<<zz<<" "<< fu[i]<<std::endl;
    }
    return fu_mem;
}

target_gradients_aniso::target_gradients_aniso(
                               af::const_ref<double> const& fo,
                               af::const_ref< std::complex<double> > const& fc,
                               af::const_ref< std::complex<double> > const& fm,
                               sym_mat3<double> const& b_cart,
                               double const& ksol,
                               double const& bsol,
                               af::const_ref<cctbx::miller::index<> > const& hkl,
                               cctbx::uctbx::unit_cell const& uc,
                               bool const& calc_grad_u,
                               bool const& calc_grad_ksol,
                               bool const& calc_grad_bsol)
{
    af::shared<double> s_mem = uc.stol_sq(hkl);
    af::const_ref<double> s = s_mem.const_ref();
    mat3<double> a = uc.fractionalization_matrix();
    sym_mat3<double> u_star = b_cart.tensor_transform(a);
    MMTBX_ASSERT(hkl.size()==fo.size());
    MMTBX_ASSERT(fo.size()==fc.size()&&fc.size()==fm.size());
    std::vector<double> f_b(fo.size());
    std::vector<double> fmodel_complex_abs(fo.size());
    af::shared<std::complex<double> > fmodel_complex(fo.size());
    double num=0.0;
    double denum=0.0;
    for(std::size_t i=0; i < fo.size(); i++) {
      double fu = fu_star(u_star, hkl[i]);
      f_b[i]=std::exp(-bsol * s[i]);
      double fmodel_abs = std::abs((fc[i] + f_b[i]*ksol * fm[i]) * fu);
      fmodel_complex_abs[i]=fmodel_abs;
      fmodel_complex[i]=(fc[i] + f_b[i]*ksol * fm[i]) * fu;
      num += fo[i] * fmodel_abs;
      denum += fmodel_abs*fmodel_abs;
    }
    //af::const_ref<std::complex<double> > fmodel_complex_ = fmodel_complex.const_ref();
    //af::shared<double> weights(fo.size(), 1.0);
    //af::const_ref<double> weights_ = weights.const_ref();
    //cctbx::xray::targets::ls_target_with_scale_k1 tgx_manager(
    //                       fo,
    //                       weights_,
    //                       fmodel_complex_,
    //                       true,
    //                       false,
    //                       0.0);
    //af::shared<std::complex<double> > d_target_d_fmodel = tgx_manager.derivatives();

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
// Matamatica proven equivalent of gxu
        //int h = mi[0];
        //int k = mi[1];
        //int l = mi[2];
        //
        //double p0 = a[0]*a[0]*h*h +   a[3]*a[3]*k*k + 2*a[3]*a[6]*k*l + a[6]*a[6]*l*l + 2*a[0]*h*(a[3]*k + a[6]*l);
        //double p1 = a[1]*a[1]*h*h +   a[4]*a[4]*k*k + 2*a[4]*a[7]*k*l + a[7]*a[7]*l*l + 2*a[1]*h*(a[4]*k + a[7]*l);
        //double p2 = a[2]*a[2]*h*h +    2*a[2]*a[5]*h*k + a[5]*a[5]*k*k + 2*(a[2]*a[8]*h + a[5]*a[8]*k)*l + a[8]*a[8]*l*l;
        //double p3 = 2*a[3]*a[4]*k*k + 2*(a[4]*a[6] + a[3]*a[7])*k*l + 2*a[6]*a[7]*l*l + 2*a[1]*h*(a[3]*k + a[6]*l) + 2*a[0]*h*(a[1]*h + a[4]*k + a[7]*l);
        //double p4 = 2*a[2]*a[3]*h*k + 2*a[3]*a[5]*k*k + 2*(a[2]*a[6]*h + (a[5]*a[6] + a[3]*a[8])*k)*l + 2*a[6]*a[8]*l*l + 2*a[0]*h*(a[2]*h + a[5]*k + a[8]*l);
        //double p5 = 2*a[2]*a[4]*h*k + 2*a[4]*a[5]*k*k + 2*(a[2]*a[7]*h + (a[5]*a[7] + a[4]*a[8])*k)*l + 2*a[7]*a[8]*l*l + 2*a[1]*h*(a[2]*h + a[5]*k + a[8]*l);

        //std::cout<< " " <<std::endl;
        //std::cout<< p0<<" =0= "<<t1*t1     <<std::endl;
        //std::cout<< p1<<" =1= "<<t2*t2     <<std::endl;
        //std::cout<< p2<<" =2= "<<t3*t3     <<std::endl;
        //std::cout<< p3<<" =3= "<<2.*t1*t2 <<std::endl;
        //std::cout<< p4<<" =4= "<<2.*t1*t3 <<std::endl;
        //std::cout<< p5<<" =5= "<<2.*t2*t3 <<std::endl;

        //gxu[0] += coeff * p0;
        //gxu[1] += coeff * p1;
        //gxu[2] += coeff * p2;
        //gxu[3] += coeff * p3;
        //gxu[4] += coeff * p4;
        //gxu[5] += coeff * p5;

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
   //tgx = tgx_manager.target();
   //std::cout<<tgx_manager.target()<<" "<< tgx<<std::endl;
   gtgx_ksol /= scale_tgx;
   gtgx_bsol /= scale_tgx;
   for(std::size_t j=0; j < 6; j++) {
     gxu[j] /= scale_tgx;
   }
}



double r_factor_aniso_fast(af::const_ref<double> const& fo,
                           af::const_ref< std::complex<double> > const& fc,
                           af::const_ref< std::complex<double> > const& fm,
                           sym_mat3<double> const& b_cart,
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
    sym_mat3<double> u_star = sym_mat3<double> (b_cart).tensor_transform(a);
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
                               sym_mat3<double> const& b_cart,
                               double const& ksol,
                               double const& bsol,
                               af::const_ref<cctbx::miller::index<> > const& hkl,
                               cctbx::uctbx::unit_cell const& uc,
                               cctbx::sgtbx::space_group const& sg,
                               af::const_ref<bool> const& gradient_flags,
                               af::const_ref<double> const& alpha,
                               af::const_ref<double> const& beta,
                               double k)
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
    sym_mat3<double> u_star = sym_mat3<double> (b_cart).tensor_transform(a);
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
}

}} // namespace mmtbx::bulk_solvent
