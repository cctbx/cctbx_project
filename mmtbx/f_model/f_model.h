#ifndef MMTBX_F_MODEL_H
#define MMTBX_F_MODEL_H

#include <mmtbx/error.h>
#include <cctbx/miller.h>
#include <cctbx/uctbx.h>
#include <cctbx/xray/targets.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/sym_mat3.h>

using namespace std;
namespace mmtbx { namespace f_model {
namespace af=scitbx::af;
using scitbx::mat3;
using scitbx::sym_mat3;

template <typename FloatType>
FloatType f_aniso_one_h(cctbx::miller::index<> const& h,
                        scitbx::sym_mat3<FloatType> u_star)
{
  FloatType arg = -2.0*scitbx::constants::pi*scitbx::constants::pi *
       (u_star[0]*h[0]*h[0] +
        u_star[1]*h[1]*h[1] +
        u_star[2]*h[2]*h[2] +
     2.*u_star[3]*h[0]*h[1] +
     2.*u_star[4]*h[0]*h[2] +
     2.*u_star[5]*h[1]*h[2]);
  if(arg > 40.0) arg=40.0; // to avoid overflow problem
  return std::exp(arg);
}

template <typename FloatType>
FloatType f_b_exp_one_h(FloatType const& ss, FloatType const& b)
{
  FloatType arg = -ss * b;
  if(arg > 40.0) arg=40.0; // to avoid overflow problem
  return std::exp(arg);
}

template <typename FloatType, typename ComplexType>
ComplexType f_model_one_h(ComplexType const& f_calc,
                          ComplexType const& f_mask,
                          FloatType const& k_sol,
                          FloatType const& b_sol,
                          FloatType const& ss,
                          cctbx::miller::index<> const& h,
                          scitbx::sym_mat3<FloatType> u_star)
{
  return
    f_aniso_one_h(h,u_star)*(f_calc+f_mask*k_sol*f_b_exp_one_h(ss,b_sol));
}

template <typename FloatType, typename ComplexType>
af::shared<ComplexType> f_model(
  af::const_ref<ComplexType> const& f_calc,
  af::const_ref<ComplexType> const& f_mask,
  sym_mat3<FloatType> const& u_star,
  FloatType const& k_sol,
  FloatType const& b_sol,
  af::const_ref<FloatType> const& ss,
  af::const_ref<cctbx::miller::index<> > const& hkl)
{
  MMTBX_ASSERT(f_calc.size()==f_mask.size());
  MMTBX_ASSERT(f_calc.size()==ss.size());
  MMTBX_ASSERT(f_calc.size()==hkl.size());
  af::shared<std::complex<FloatType> > result(f_calc.size());
  for(std::size_t i=0; i < f_calc.size(); i++) {
    result[i] = f_model_one_h(f_calc[i], f_mask[i], k_sol, b_sol, ss[i], hkl[i],
      u_star);
  }
  return result;
}

template <typename FloatType=double,
          typename ComplexType=std::complex<double> >
class core
{
  public:
    af::shared<ComplexType>                f_calc;
    af::shared<ComplexType>                f_mask;
    FloatType                              k_sol;
    FloatType                              b_sol;
    af::shared<ComplexType>                f_part_base;
    af::shared<ComplexType>                f_part;
    FloatType                              k_part;
    FloatType                              b_part;
    scitbx::sym_mat3<FloatType>            u_star;
    af::shared<cctbx::miller::index<> >    hkl;
    cctbx::uctbx::unit_cell                uc;
    af::shared<FloatType>                  ss, f_aniso, f_b_sol;
    af::shared<ComplexType>                f_model, f_bulk;
    mat3<FloatType>                        a;

    core() {}

    core(af::shared<ComplexType>                const& f_calc_,
         af::shared<ComplexType>                const& f_mask_,
         FloatType                              const& k_sol_,
         FloatType                              const& b_sol_,
         af::shared<ComplexType>                const& f_part_base_,
         FloatType                              const& k_part_,
         FloatType                              const& b_part_,
         scitbx::sym_mat3<FloatType>            const& u_star_,
         af::shared<cctbx::miller::index<> >    const& hkl_,
         cctbx::uctbx::unit_cell                const& uc_,
         af::shared<FloatType>                  const& ss_)
    :
      f_calc(f_calc_), f_mask(f_mask_), k_sol(k_sol_),
      b_sol(b_sol_), u_star(u_star_), hkl(hkl_),uc(uc_), ss(ss_),
      f_aniso(hkl_.size(), af::init_functor_null<FloatType>()),
      f_bulk(hkl_.size(), af::init_functor_null<ComplexType>()),
      f_model(hkl_.size(), af::init_functor_null<ComplexType>()),
      f_b_sol(hkl_.size(), af::init_functor_null<FloatType>()),
      a(uc_.fractionalization_matrix()), k_part(k_part_), b_part(b_part_),
      f_part_base(f_part_base_),
      f_part(f_part_base_.size(), af::init_functor_null<ComplexType>())
    {
      MMTBX_ASSERT(f_calc.size() == f_mask.size());
      MMTBX_ASSERT(f_calc.size() == hkl.size()   );
      MMTBX_ASSERT(f_calc.size() == f_calc.size());
      MMTBX_ASSERT(f_calc.size() == f_calc.size());
      MMTBX_ASSERT(f_calc.size() == f_part_base_.size());
      FloatType*   f_aniso_  = f_aniso.begin();
      FloatType*   f_b_sol_  = f_b_sol.begin();
      FloatType*   ss__      = ss.begin();
      ComplexType* f_bulk_   = f_bulk.begin();
      ComplexType* f_model_  = f_model.begin();
      ComplexType* f_calc__  = f_calc.begin();
      ComplexType* f_mask__  = f_mask.begin();
      ComplexType* f_part_   = f_part.begin();
      ComplexType* f_part_base__ = f_part_base.begin();
      for(std::size_t i=0; i < hkl.size(); i++) {
        f_aniso_[i] = f_aniso_one_h(hkl[i], u_star);
        f_b_sol_[i] = f_b_exp_one_h(ss__[i], b_sol);
        f_bulk_[i]  = f_k_bexp_f_one_h(f_mask__[i], f_b_sol_[i], k_sol);
        f_part_[i]  = f_k_bexp_f_one_h(f_part_base__[i],
                                       f_b_exp_one_h(ss__[i], b_part), k_part);
        f_model_[i] = f_model_one_h(f_calc__[i], f_bulk_[i], f_part_[i], f_aniso_[i]);
      }
    }

    protected:
       ComplexType f_k_bexp_f_one_h(ComplexType const& f_h,
                                    FloatType   const& f_bexp_h,
                                    FloatType   const& k)
       {
         return f_h * f_bexp_h * k;
       }
       ComplexType f_model_one_h(ComplexType const& f_calc_h,
                                 ComplexType const& f_bulk_h,
                                 ComplexType const& f_part_h,
                                 FloatType   const& f_aniso_h)
       {
         return f_aniso_h * (f_calc_h + f_bulk_h + f_part_h);
       }
};

}} // namespace mmtbx::f_model

#endif // MMTBX_F_MODEL_H
