#include <mmtbx/bulk_solvent/bulk_solvent.h>
#include <cctbx/xray/targets.h>

namespace mmtbx { namespace bulk_solvent {


double fu_star(sym_mat3<double> const& u_star,
               cctbx::miller::index<> const& mi)
{
    double arg = -0.25 * (u_star[0]*mi[0]*mi[0] +
                          u_star[1]*mi[1]*mi[1] +
                          u_star[2]*mi[2]*mi[2] +
                       2.*u_star[3]*mi[0]*mi[1] +
                       2.*u_star[4]*mi[0]*mi[2] +
                       2.*u_star[5]*mi[1]*mi[2]);
    if(arg > 40.0) arg=40.0; // to avoid overflow problem
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
      fu[i] = fu_star(u_star, hkl[i]);
    }
    return fu_mem;
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
      bool cs = cf[i];
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
