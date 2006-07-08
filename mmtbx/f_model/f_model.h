#ifndef MMTBX_F_MODEL_H
#define MMTBX_F_MODEL_H

#include <scitbx/array_family/shared.h>
#include <mmtbx/error.h>
#include <cctbx/miller.h>
#include <cctbx/xray/targets.h>

using namespace std;
namespace mmtbx { namespace f_model {
namespace af=scitbx::af;
using scitbx::mat3;
using scitbx::sym_mat3;

template <typename FloatType=double,
          typename ComplexType=std::complex<double> >
class core
{
  public:
    af::shared<ComplexType>                f_calc;
    af::shared<ComplexType>                f_mask;
    scitbx::sym_mat3<FloatType>            b_cart;
    FloatType                              k_sol;
    FloatType                              b_sol;
    af::const_ref<cctbx::miller::index<> > hkl;
    cctbx::uctbx::unit_cell                uc;
    af::shared<FloatType>                  ss, fb_cart;
    af::shared<ComplexType>                f_model, f_bulk;
    mat3<FloatType> a;

    core() {}

    core(af::shared<ComplexType>                const& f_calc_,
         af::shared<ComplexType>                const& f_mask_,
         scitbx::sym_mat3<FloatType>            const& b_cart_,
         FloatType                              const& k_sol_,
         FloatType                              const& b_sol_,
         af::const_ref<cctbx::miller::index<> > const& hkl_,
         cctbx::uctbx::unit_cell                const& uc_,
         af::shared<FloatType>                  const& ss_)
    :
      f_calc(f_calc_), f_mask(f_mask_),   k_sol(k_sol_),
      b_sol(b_sol_), b_cart(b_cart_), hkl(hkl_),uc(uc_), ss(ss_),
      fb_cart(hkl_.size(), af::init_functor_null<FloatType>()),
      f_bulk(hkl_.size(), af::init_functor_null<ComplexType>()),
      f_model(hkl_.size(), af::init_functor_null<ComplexType>()),
      a(uc_.fractionalization_matrix())
    {
      MMTBX_ASSERT(f_calc.size() == f_mask.size());
      MMTBX_ASSERT(f_calc.size() == hkl.size()   );
      MMTBX_ASSERT(f_calc.size() == f_calc.size());
      MMTBX_ASSERT(f_calc.size() == f_calc.size());
      FloatType*   fb_cart_  = fb_cart.begin();
      FloatType*   ss__       = ss.begin();
      ComplexType* f_bulk_   = f_bulk.begin();
      ComplexType* f_model_  = f_model.begin();
      ComplexType* f_calc__  = f_calc.begin();
      ComplexType* f_mask__  = f_mask.begin();
      for(std::size_t i=0; i < hkl.size(); i++) {
          fb_cart_[i] = fb_cart_one_h(b_cart, hkl[i]);
          f_bulk_[i]  = f_bulk_one_h(f_mask__[i], ss__[i]);
          f_model_[i] = f_model_one_h(f_calc__[i], f_bulk_[i], fb_cart_[i]);
      }
    }

    protected:
       FloatType fb_cart_one_h(sym_mat3<FloatType> const& b_cart,
                               cctbx::miller::index<> const& h)
       {
           sym_mat3<double> u_star =
                                 sym_mat3<double> (b_cart).tensor_transform(a);
           FloatType arg = -0.25 * (u_star[0]*h[0]*h[0] +
                                    u_star[1]*h[1]*h[1] +
                                    u_star[2]*h[2]*h[2] +
                                 2.*u_star[3]*h[0]*h[1] +
                                 2.*u_star[4]*h[0]*h[2] +
                                 2.*u_star[5]*h[1]*h[2]);
           if(arg > std::log(5.0)) return 1.0;
           return std::exp(arg);
       }
       ComplexType f_bulk_one_h(ComplexType const& f_mask_h,
                                FloatType   const& ss_h)
       {
           return f_mask_h * std::exp(-ss_h * b_sol) * k_sol;
       }
       ComplexType f_model_one_h(ComplexType const& f_calc_h,
                                 ComplexType const& f_bulk_h,
                                 FloatType   const& fb_cart_h)
       {
           return fb_cart_h * (f_calc_h + f_bulk_h);
       }
};

template <typename FloatType=double,
          typename ComplexType=std::complex<double> >
class d_fmodel_d_kb_sol
{
  public:
    af::shared<ComplexType> d_fmodel_d_ksol, d_fmodel_d_bsol;

    d_fmodel_d_kb_sol() {}

    d_fmodel_d_kb_sol(core<FloatType, ComplexType> core_)
    :
      d_fmodel_d_ksol(core_.ss.size(), af::init_functor_null<ComplexType>()),
      d_fmodel_d_bsol(core_.ss.size(), af::init_functor_null<ComplexType>())
    {
      FloatType    k_sol     = core_.k_sol;
      FloatType    b_sol     = core_.b_sol;
      FloatType*   ss        = core_.ss.begin();
      ComplexType* f_calc    = core_.f_calc.begin();
      ComplexType* f_mask    = core_.f_mask.begin();
      ComplexType* f_model   = core_.f_model.begin();
      ComplexType* d_fmodel_d_ksol_ = d_fmodel_d_ksol.begin();
      ComplexType* d_fmodel_d_bsol_ = d_fmodel_d_bsol.begin();
      af::const_ref<cctbx::miller::index<> > hkl = core_.hkl;
      for(std::size_t i=0; i < hkl.size(); i++) {
          FloatType uvs_plus_usv = std::real(f_calc[i]*std::conj(f_mask[i])+
                                   f_mask[i]*std::conj(f_calc[i]));
          FloatType mod_v_sq = std::abs(f_mask[i]) * std::abs(f_mask[i]);
          FloatType f_b = std::exp(-b_sol * ss[i]);
          FloatType coeff =
                  (uvs_plus_usv+f_b*k_sol*mod_v_sq*2)/(std::abs(f_model[i])*2);
          d_fmodel_d_ksol_[i] = coeff * f_b;
          d_fmodel_d_bsol_[i] = - coeff * f_b * k_sol * ss[i];
      }
    }
};

template <typename FloatType=double,
          typename ComplexType=std::complex<double> >
class ls_target_and_kbu_gradients
{
  public:
    FloatType d_target_d_ksol, d_target_d_bsol, target;
    scitbx::sym_mat3<FloatType> d_target_d_u;
    af::shared<FloatType> f_obs;
    bool calc_grad_u, calc_grad_ksol, calc_grad_bsol;
    core<FloatType, ComplexType> core_data;

    ls_target_and_kbu_gradients() {}

    ls_target_and_kbu_gradients(core<FloatType, ComplexType> const& core_,
                                af::shared<FloatType>        const& f_obs_,
                                bool const& calc_grad_u_,
                                bool const& calc_grad_ksol_,
                                bool const& calc_grad_bsol_)
    :
      d_target_d_ksol(0.0), d_target_d_bsol(0.0), d_target_d_u(0.0),
      f_obs(f_obs_),calc_grad_u(calc_grad_u_), calc_grad_ksol(calc_grad_ksol_),
      calc_grad_bsol(calc_grad_bsol_),core_data(core_)
    {
      af::shared<FloatType> weights(f_obs_.size(), 1.0);
      cctbx::xray::targets::ls_target_with_scale_k1 tgx_manager(
                                                     f_obs.const_ref(),
                                                     weights.const_ref(),
                                                     core_.f_model.const_ref(),
                                                     true,
                                                     false,
                                                     0.0);
      target = tgx_manager.target();
      af::shared<ComplexType> d_target_d_fmodel = tgx_manager.derivatives();
      d_fmodel_d_kb_sol<FloatType, ComplexType> d_fmodel_d_kb_sol_manager(core_data);
      ComplexType* f_model   = core_data.f_model.begin();
      ComplexType* d_target_d_fmodel_ = d_target_d_fmodel.begin();
      ComplexType* d_fmodel_d_ksol_ =
                             d_fmodel_d_kb_sol_manager.d_fmodel_d_ksol.begin();
      ComplexType* d_fmodel_d_bsol_ =
                             d_fmodel_d_kb_sol_manager.d_fmodel_d_bsol.begin();
      for(std::size_t i=0; i < f_obs.size(); i++) {
          ComplexType correcting_factor = f_model[i] / std::conj(f_model[i]);
          if(calc_grad_ksol) {
             ComplexType result =
               d_target_d_fmodel_[i] * d_fmodel_d_ksol_[i] * correcting_factor;
             d_target_d_ksol += result.real();
             MMTBX_ASSERT(std::abs(result.imag()) < 1.e-6);
          }
          if(calc_grad_bsol) {
             ComplexType result =
               d_target_d_fmodel_[i] * d_fmodel_d_bsol_[i] * correcting_factor;
             d_target_d_bsol += result.real();
             MMTBX_ASSERT(std::abs(result.imag()) < 1.e-6);
          }
      }
    }
};

}} // namespace mmtbx::f_model

#endif // MMTBX_F_MODEL_H
