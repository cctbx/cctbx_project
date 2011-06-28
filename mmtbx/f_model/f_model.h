#ifndef MMTBX_F_MODEL_H
#define MMTBX_F_MODEL_H

#include <mmtbx/error.h>
#include <cctbx/adptbx.h>
#include <cctbx/xray/targets.h>

namespace mmtbx { namespace f_model {
using namespace std;
namespace af=scitbx::af;
using scitbx::mat3;
using scitbx::sym_mat3;

const int max_n_shells = 10;

template <typename FloatType>
FloatType f_aniso_one_h(cctbx::miller::index<> const& h,
                        scitbx::sym_mat3<FloatType> u_star)
{
  return cctbx::adptbx::debye_waller_factor_u_star(
    h, u_star,
    /*exp_arg_limit*/ 40., /*truncate_exp_arg*/ true);
}

template <typename FloatType>
FloatType f_b_exp_one_h(FloatType const& ss, FloatType const& b)
{
  return cctbx::adptbx::debye_waller_factor_b_iso(
    /*stol_sq*/ ss, /*b_iso*/ b,
    /*exp_arg_limit*/ 40., /*truncate_exp_arg*/ true);
}

template <typename FloatType=double,
          typename ComplexType=std::complex<double> >
class core
{
  public:
    af::shared<ComplexType>                f_calc;
    FloatType                              b_sol;
    af::shared<ComplexType>                f_part1;
    af::shared<ComplexType>                f_part2;
    scitbx::sym_mat3<FloatType>            u_star;
    af::shared<cctbx::miller::index<> >    hkl;
    cctbx::uctbx::unit_cell                uc;
    af::shared<FloatType>                  ss, f_aniso, f_b_sol;
    af::shared<ComplexType>                f_model, f_bulk;
    mat3<FloatType>                        a;

    int n_shells() const
    {
      MMTBX_ASSERT(k_sols_.size()==shell_f_masks_.size());
      return k_sols_.size();
    }

    FloatType k_sol(unsigned short i) const
    {
      return k_sols_[i];
    }

    af::small<FloatType,max_n_shells> k_sols() const
    {
      return k_sols_;
    }

    af::shared<ComplexType> f_mask() const
    {
      MMTBX_ASSERT(shell_f_masks_.size()==1U);
      return shell_f_masks_[0];
    }

    af::shared<ComplexType>  shell_f_mask(unsigned short i) const
    {
      return shell_f_masks_[i];
    }

    af::small<af::shared<ComplexType>,max_n_shells> shell_f_masks() const
    {
      return shell_f_masks_;
    }


    core() {}

    core(af::shared<ComplexType>                const& f_calc_,
         af::small< af::shared<ComplexType>, max_n_shells> const& f_masks,
         af::small<FloatType, max_n_shells>     const& k_sols,
         FloatType                              const& b_sol_,
         af::shared<ComplexType>                const& f_part1_,
         af::shared<ComplexType>                const& f_part2_,
         scitbx::sym_mat3<FloatType>            const& u_star_,
         af::shared<cctbx::miller::index<> >    const& hkl_,
         cctbx::uctbx::unit_cell                const& uc_,
         af::shared<FloatType>                  const& ss_)
    :
      f_calc(f_calc_),
      b_sol(b_sol_), u_star(u_star_), hkl(hkl_),uc(uc_), ss(ss_),
      f_aniso(hkl_.size(), af::init_functor_null<FloatType>()),
      f_bulk(hkl_.size(), af::init_functor_null<ComplexType>()),
      f_model(hkl_.size(), af::init_functor_null<ComplexType>()),
      f_b_sol(hkl_.size(), af::init_functor_null<FloatType>()),
      a(uc_.fractionalization_matrix()),
      f_part1(f_part1_), f_part2(f_part2_)
    {
      MMTBX_ASSERT(f_calc.size() == hkl.size());
      for(short j=0; j<f_masks.size(); ++j)
        MMTBX_ASSERT(f_calc.size() == f_masks[j].size());
      MMTBX_ASSERT(f_calc.size() == f_part1_.size());
      MMTBX_ASSERT(f_calc.size() == f_part2_.size());
      MMTBX_ASSERT(k_sols.size() == f_masks.size());
      k_sols_ = k_sols;
      shell_f_masks_ = f_masks;
      MMTBX_ASSERT( this->n_shells() >= 1U );
      for(short j=0; j<shell_f_masks_.size(); ++j)
        MMTBX_ASSERT(f_calc.size() == this->shell_f_masks_[j].size());
      FloatType*   f_aniso_  = f_aniso.begin();
      FloatType*   f_b_sol_  = f_b_sol.begin();
      FloatType*   ss__      = ss.begin();
      ComplexType* f_bulk_   = f_bulk.begin();
      ComplexType* f_model_  = f_model.begin();
      ComplexType* f_calc__  = f_calc.begin();
      ComplexType* f_part1__ = f_part1.begin();
      ComplexType* f_part2__ = f_part2.begin();
      scitbx::af::small<ComplexType*,max_n_shells> f_masks__;
      for(short j=0; j<shell_f_masks_.size(); ++j)
        f_masks__.push_back(shell_f_masks_[j].begin());
      for(std::size_t i=0; i < hkl.size(); i++) {
        f_aniso_[i] = f_aniso_one_h(hkl[i], u_star);
        f_b_sol_[i] = f_b_exp_one_h(ss__[i], b_sol);
        // TODO: optimize
        f_bulk_[i] = f_k_bexp_f_one_h(f_masks__[0][i], f_b_sol_[i], k_sols_[0]);
        for(short j=1; j<k_sols_.size(); ++j)
          f_bulk_[i] +=f_k_bexp_f_one_h(f_masks__[j][i],f_b_sol_[i],k_sols_[j]);
        f_model_[i] = f_model_one_h(
          f_calc__[i], f_bulk_[i], f_part1_[i], f_part2__[i], f_aniso_[i]);
      }
    }

    protected:
       ComplexType f_k_bexp_f_one_h(ComplexType const& f_h,
                                    FloatType   const& f_bexp_h,
                                    FloatType   const& k)
       {
         return k * f_bexp_h * f_h;
       }
       ComplexType f_model_one_h(ComplexType const& f_calc_h,
                                 ComplexType const& f_bulk_h,
                                 ComplexType const& f_part1_h,
                                 ComplexType const& f_part2_h,
                                 FloatType   const& f_aniso_h)
       {
         return f_aniso_h * (f_calc_h + f_bulk_h + f_part1_h + f_part2_h);
       }

       scitbx::af::small<af::shared<ComplexType>, max_n_shells> shell_f_masks_;
       scitbx::af::small<FloatType,max_n_shells> k_sols_;
};

}} // namespace mmtbx::f_model

#endif // MMTBX_F_MODEL_H
