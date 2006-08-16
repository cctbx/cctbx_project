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
    FloatType                              overall_scale;
    af::const_ref<cctbx::miller::index<> > hkl;
    cctbx::uctbx::unit_cell                uc;
    af::shared<FloatType>                  ss, fb_cart;
    af::shared<ComplexType>                f_model, f_bulk;
    mat3<FloatType>                        a;

    core() {}

    core(af::shared<ComplexType>                const& f_calc_,
         af::shared<ComplexType>                const& f_mask_,
         scitbx::sym_mat3<FloatType>            const& b_cart_,
         FloatType                              const& k_sol_,
         FloatType                              const& b_sol_,
         FloatType                              const& overall_scale_,
         af::const_ref<cctbx::miller::index<> > const& hkl_,
         cctbx::uctbx::unit_cell                const& uc_,
         af::shared<FloatType>                  const& ss_)
    :
      f_calc(f_calc_), f_mask(f_mask_), k_sol(k_sol_),
      b_sol(b_sol_), b_cart(b_cart_), hkl(hkl_),uc(uc_), ss(ss_),
      fb_cart(hkl_.size(), af::init_functor_null<FloatType>()),
      f_bulk(hkl_.size(), af::init_functor_null<ComplexType>()),
      f_model(hkl_.size(), af::init_functor_null<ComplexType>()),
      a(uc_.fractionalization_matrix()),overall_scale(overall_scale_)
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
          fb_cart_[i] = fb_cart_one_h(hkl[i]);
          f_bulk_[i]  = f_bulk_one_h(f_mask__[i], ss__[i]);
          f_model_[i] = f_model_one_h(f_calc__[i], f_bulk_[i], fb_cart_[i]);
      }
    }

    protected:
       FloatType fb_cart_one_h(cctbx::miller::index<> const& h)
       {
           sym_mat3<double> u_star =
                                 sym_mat3<double> (b_cart).tensor_transform(a);
           FloatType arg = -0.25 * (u_star[0]*h[0]*h[0] +
                                    u_star[1]*h[1]*h[1] +
                                    u_star[2]*h[2]*h[2] +
                                 2.*u_star[3]*h[0]*h[1] +
                                 2.*u_star[4]*h[0]*h[2] +
                                 2.*u_star[5]*h[1]*h[2]);
           if(arg > std::log(5.0)) return 1.0; // to avoid overflow problem
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
           return overall_scale * fb_cart_h * (f_calc_h + f_bulk_h);
       }
};

}} // namespace mmtbx::f_model

#endif // MMTBX_F_MODEL_H
