#ifndef MMTBX_BULK_SOLVENT_BULK_SOLVENT_H
#define MMTBX_BULK_SOLVENT_BULK_SOLVENT_H

#include <mmtbx/error.h>
#include <mmtbx/import_scitbx_af.h>
#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <mmtbx/f_model/f_model.h>

namespace mmtbx { namespace detail {

template <typename FloatType=double,
          typename ComplexType=std::complex<double> >
class d_f_model_d_k_sol_and_d_b_sol_one_h
{
public:
  d_f_model_d_k_sol_and_d_b_sol_one_h(
    f_model::core<FloatType,ComplexType> const& fm, std::size_t i)
  {
    grad_b_sol = 0;
    grad_k_sols.resize(fm.n_shells(), 0.);
    ComplexType f_mdl = fm.f_model[i];
    FloatType f_model_abs = std::abs(f_mdl);
    if(f_model_abs > 0){
      FloatType f_b = fm.f_b_sol[i];
      ComplexType k_f_mask(0.,0.);
      for(unsigned short j=0; j<fm.n_shells(); ++j)
      {
        ComplexType f_m = fm.shell_f_mask(j)[i];
        FloatType ksol = fm.k_sol(j);
        FloatType uvs_plus_usv = std::real(f_mdl*std::conj(f_m)+f_m*std::conj(f_mdl));
        k_f_mask += ksol * f_m;
        FloatType theta = (uvs_plus_usv)/(f_model_abs*2);
        FloatType coeff = theta * fm.f_aniso[i];
        grad_k_sols[j] =  coeff * f_b;
      }
      FloatType uvs_plus_usv_bsol =
        std::real(f_mdl*std::conj(k_f_mask)+k_f_mask*std::conj(f_mdl));
      grad_b_sol =-uvs_plus_usv_bsol*fm.f_aniso[i]*f_b*fm.ss[i]/(f_model_abs*2);
    }
  }
  FloatType grad_b_sol;
  af::small<FloatType,mmtbx::f_model::max_n_shells> grad_k_sols;
};

template <typename FloatType, typename ComplexType>
scitbx::af::tiny<FloatType, 6>
d_f_model_d_u_star_one_h(f_model::core<FloatType,ComplexType> const& fm,
                         std::size_t i)
{
  scitbx::af::tiny<FloatType, 6> result;
  FloatType f_model_abs = std::abs(fm.f_model[i]);
  FloatType minus_two_pi = -2.0*scitbx::constants::pi*scitbx::constants::pi;
  FloatType coeff = f_model_abs * minus_two_pi;
  cctbx::miller::index<> const& mi = fm.hkl[i];
  result = scitbx::af::tiny<FloatType, 6> (
    coeff * mi[0]*mi[0],
    coeff * mi[1]*mi[1],
    coeff * mi[2]*mi[2],
    coeff * 2.*mi[0]*mi[1],
    coeff * 2.*mi[0]*mi[2],
    coeff * 2.*mi[1]*mi[2]);
  return result;
};

template <typename FloatType=double,
          typename ComplexType=std::complex<double> >
class one_h_ls
{
public:
  one_h_ls(FloatType const& fo,
           f_model::core<FloatType,ComplexType> const& fm,
           std::size_t i,
           FloatType const& scale,
           bool const& compute_k_sol_grad,
           bool const& compute_b_sol_grad,
           bool const& compute_u_star_grad)
  {
    grad_k_sols.resize(fm.n_shells(), 0.);
    FloatType f_model_abs = std::abs(fm.f_model[i]);
    diff = fo - scale * f_model_abs;
    FloatType mtsd = 0;
    if(compute_k_sol_grad || compute_b_sol_grad || compute_u_star_grad) {
      mtsd = -2. * scale * diff;
    }
    if(compute_k_sol_grad || compute_b_sol_grad) {
      d_f_model_d_k_sol_and_d_b_sol_one_h<FloatType,ComplexType> kbsol_grads =
        d_f_model_d_k_sol_and_d_b_sol_one_h<FloatType,ComplexType> (fm,i);
      for(unsigned short j=0; j<grad_k_sols.size(); ++j)
      {
        grad_k_sols[j] = mtsd*kbsol_grads.grad_k_sols[j];
      }
      grad_b_sol = mtsd * kbsol_grads.grad_b_sol;
    }
    if(compute_u_star_grad) {
      scitbx::af::tiny<FloatType, 6> usg = d_f_model_d_u_star_one_h(fm,i);
      for(std::size_t j=0; j<6; j++) grad_u_star[j] = mtsd * usg[j];
    }
  }

  one_h_ls(FloatType const& fo,
           f_model::core<FloatType,ComplexType> const& fm1,
           f_model::core<FloatType,ComplexType> const& fm2,
           FloatType const& twin_fraction,
           std::size_t i,
           FloatType const& scale,
           bool const& compute_k_sol_grad,
           bool const& compute_b_sol_grad,
           bool const& compute_u_star_grad)
  {
    MMTBX_ASSERT( fm1.n_shells() == fm2.n_shells() );
    grad_k_sols.resize(fm1.n_shells(), 0.);
    FloatType f_model_abs1 = std::abs(fm1.f_model[i]);
    FloatType f_model_abs2 = std::abs(fm2.f_model[i]);
    FloatType f_model_abs = std::sqrt(
      (1-twin_fraction)*f_model_abs1*f_model_abs1+
         twin_fraction *f_model_abs2*f_model_abs2);
    diff = fo - scale * f_model_abs;
    FloatType mtsd = 0;
    if(compute_k_sol_grad || compute_b_sol_grad || compute_u_star_grad) {
      mtsd = -2. * scale * diff / (2 * f_model_abs);
    }
    FloatType fmi1 = 2.*std::abs(fm1.f_model[i]);
    FloatType fmi2 = 2.*std::abs(fm2.f_model[i]);
    if(compute_k_sol_grad || compute_b_sol_grad) {
      d_f_model_d_k_sol_and_d_b_sol_one_h<FloatType,ComplexType> kbsol_grads1 =
        d_f_model_d_k_sol_and_d_b_sol_one_h<FloatType,ComplexType> (fm1,i);
      d_f_model_d_k_sol_and_d_b_sol_one_h<FloatType,ComplexType> kbsol_grads2 =
        d_f_model_d_k_sol_and_d_b_sol_one_h<FloatType,ComplexType> (fm2,i);
      for(unsigned short j=0; j<grad_k_sols.size(); ++j)
      {
        grad_k_sols[j] = mtsd *
        ((1-twin_fraction)*kbsol_grads1.grad_k_sols[j]*fmi1+
            twin_fraction *kbsol_grads2.grad_k_sols[j]*fmi2);
      }
      grad_b_sol = mtsd *
        ((1-twin_fraction)*kbsol_grads1.grad_b_sol*fmi1+
            twin_fraction *kbsol_grads2.grad_b_sol*fmi2);
    }
    if(compute_u_star_grad) {
      scitbx::af::tiny<FloatType, 6> usg1 = d_f_model_d_u_star_one_h(fm1,i);
      scitbx::af::tiny<FloatType, 6> usg2 = d_f_model_d_u_star_one_h(fm2,i);
      for(std::size_t j=0; j<6; j++) {
        grad_u_star[j] = mtsd * ((1-twin_fraction)*usg1[j]*fmi1+
                                    twin_fraction *usg2[j]*fmi2);
      }
    }
  }

  FloatType diff, grad_b_sol;
  scitbx::af::tiny<FloatType, 6> grad_u_star;
  scitbx::af::small<FloatType,mmtbx::f_model::max_n_shells> grad_k_sols;
};

template <typename FloatType=double,
          typename ComplexType=std::complex<double> >
class overall_scale
{
public:
  overall_scale() {}
  overall_scale(af::const_ref<FloatType> const& fo,
                f_model::core<FloatType,ComplexType> const& fm)
  {
    FloatType num=0.;
    FloatType denum=0.;
    sum_fo_sq = 0.;
    for(std::size_t i=0; i < fo.size(); i++) {
      FloatType f_model_abs = std::abs(fm.f_model[i]);
      num += fo[i] * f_model_abs;
      denum += f_model_abs * f_model_abs;
      sum_fo_sq += fo[i]*fo[i];
    }
    scale = (denum == 0 ? 0 : num/denum);
  }
  overall_scale(af::const_ref<FloatType> const& fo,
                f_model::core<FloatType,ComplexType> const& fm1,
                f_model::core<FloatType,ComplexType> const& fm2,
                FloatType const& twin_fraction)
  {
    FloatType num=0.;
    FloatType denum=0.;
    sum_fo_sq = 0.;
    for(std::size_t i=0; i < fo.size(); i++) {
      FloatType f_model_abs1 = std::abs(fm1.f_model[i]);
      FloatType f_model_abs2 = std::abs(fm2.f_model[i]);
      FloatType f_model_abs = std::sqrt(
        (1-twin_fraction)*f_model_abs1*f_model_abs1+
           twin_fraction *f_model_abs2*f_model_abs2);
      num += fo[i] * f_model_abs;
      denum += f_model_abs * f_model_abs;
      sum_fo_sq += fo[i]*fo[i];
    }
    scale = (denum == 0 ? 0 : num/denum);
  }
  FloatType scale, sum_fo_sq;
};

}}

namespace mmtbx { namespace bulk_solvent {

using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;

template <typename FloatType=double,
          typename ComplexType=std::complex<double> >
class bulk_solvent_and_aniso_scale_target_and_grads_ls
{
public:
  bulk_solvent_and_aniso_scale_target_and_grads_ls() {}

  bulk_solvent_and_aniso_scale_target_and_grads_ls(
    f_model::core<FloatType,ComplexType> const& fm1,
    f_model::core<FloatType,ComplexType> const& fm2,
    FloatType const& twin_fraction,
    af::const_ref<FloatType> const& fo,
    bool const& compute_k_sol_grad,
    bool const& compute_b_sol_grad,
    bool const& compute_u_star_grad)
  {
    MMTBX_ASSERT(fo.size() == fm1.f_calc.size());
    MMTBX_ASSERT(fo.size() == fm2.f_calc.size());
    MMTBX_ASSERT(fm1.n_shells() == fm2.n_shells());
    grad_k_sols_.resize(fm1.n_shells(),0.);
    overall_scale =
      detail::overall_scale<FloatType,ComplexType>(fo, fm1, fm2, twin_fraction);
    grad_u_star_ = scitbx::af::tiny<FloatType, 6>(0,0,0,0,0,0);
    target_ = 0., grad_b_sol_ = 0.;
    for(std::size_t i=0; i < fo.size(); i++) {
      detail::one_h_ls<FloatType, ComplexType> one_h =
        detail::one_h_ls<FloatType,ComplexType>(fo[i],fm1,fm2,twin_fraction,i,
          overall_scale.scale, compute_k_sol_grad,compute_b_sol_grad,
          compute_u_star_grad);
      target_ += one_h.diff * one_h.diff;
      if(compute_u_star_grad) {
        for(std::size_t j=0; j<6; j++) grad_u_star_[j] += one_h.grad_u_star[j];
      }
      if(compute_k_sol_grad || compute_b_sol_grad) {
        for(unsigned short j=0; j<grad_k_sols_.size(); ++j)
        {
          grad_k_sols_[j] += one_h.grad_k_sols[j];
        }
        grad_b_sol_ += one_h.grad_b_sol;
      }
    }
    for(std::size_t i=0; i < 6; i++) grad_u_star_[i] /= overall_scale.sum_fo_sq;
  }

  bulk_solvent_and_aniso_scale_target_and_grads_ls(
    f_model::core<FloatType,ComplexType> const& fm,
    af::const_ref<FloatType> const& fo,
    bool const& compute_k_sol_grad,
    bool const& compute_b_sol_grad,
    bool const& compute_u_star_grad)
  {
    MMTBX_ASSERT(fo.size() == fm.f_calc.size());
    grad_k_sols_.resize(fm.n_shells(),0.);
    overall_scale = detail::overall_scale<FloatType,ComplexType>(fo, fm);
    grad_u_star_ = scitbx::af::tiny<FloatType, 6>(0,0,0,0,0,0);
    target_ = 0., grad_b_sol_ = 0.;
    for(std::size_t i=0; i < fo.size(); i++) {
      detail::one_h_ls<FloatType, ComplexType> one_h =
        detail::one_h_ls<FloatType,ComplexType>(fo[i],fm,i,overall_scale.scale,
          compute_k_sol_grad,compute_b_sol_grad,compute_u_star_grad);
      target_ += one_h.diff * one_h.diff;
      if(compute_u_star_grad) {
        for(std::size_t j=0; j<6; j++) grad_u_star_[j] += one_h.grad_u_star[j];
      }
      if(compute_k_sol_grad || compute_b_sol_grad) {
        for(unsigned short j=0; j<grad_k_sols_.size(); ++j)
        {
          grad_k_sols_[j] += one_h.grad_k_sols[j];
        }
        grad_b_sol_ += one_h.grad_b_sol;
      }
    }
    for(std::size_t i=0; i < 6; i++) grad_u_star_[i] /= overall_scale.sum_fo_sq;
  }

  FloatType target() { return target_/overall_scale.sum_fo_sq; }
  scitbx::af::tiny<FloatType, 6> grad_u_star() { return grad_u_star_; }

  FloatType grad_k_sol()
  {
    MMTBX_ASSERT( grad_k_sols_.size()==1U );
    return grad_k_sols_[0]/overall_scale.sum_fo_sq;
  }

  FloatType grad_k_sol(unsigned short j) const
  {
    return grad_k_sols_[j]/overall_scale.sum_fo_sq;
  }

  scitbx::af::small<FloatType,mmtbx::f_model::max_n_shells> grad_k_sols()
    const
  {
    scitbx::af::small<FloatType,mmtbx::f_model::max_n_shells>
      result(grad_k_sols_.size(),0.);
    for(unsigned short j=0; j<grad_k_sols_.size(); ++j)
    {
      result[j] = grad_k_sols_[j]/overall_scale.sum_fo_sq;
    }
    return result;
  }

  FloatType grad_b_sol() { return grad_b_sol_/overall_scale.sum_fo_sq; }

private:
  FloatType target_, grad_b_sol_;
  scitbx::af::tiny<FloatType, 6> grad_u_star_;
  detail::overall_scale<FloatType, ComplexType> overall_scale;
  scitbx::af::small<FloatType,mmtbx::f_model::max_n_shells> grad_k_sols_;
};

//------------------------------------------------------------------------------
template <typename FloatType, typename ComplexType>
FloatType
scale(af::const_ref<FloatType> const& fo,
      af::const_ref<ComplexType> const& fc)
{
    MMTBX_ASSERT(fo.size()==fc.size());
    FloatType num=0.0;
    FloatType denum=0.0;
    for(std::size_t i=0; i < fo.size(); i++) {
      FloatType fc_abs = std::abs(fc[i]);
      num += fo[i] * fc_abs;
      denum += fc_abs * fc_abs;
    }
    return (denum == 0 ? 0 : num/denum);
};

template <typename FloatType>
FloatType
scale(af::const_ref<FloatType> const& fo,
      af::const_ref<FloatType> const& fc)
{
    MMTBX_ASSERT(fo.size()==fc.size());
    FloatType num=0.0;
    FloatType denum=0.0;
    for(std::size_t i=0; i < fo.size(); i++) {
      num += fo[i] * fc[i];
      denum += fc[i] * fc[i];
    }
    return (denum == 0 ? 0 : num/denum);
};

template <typename FloatType, typename ComplexType>
FloatType
scale(
  af::const_ref<FloatType> const& fo,
  af::const_ref< std::complex<ComplexType> > const& fc1,
  af::const_ref< std::complex<ComplexType> > const& fc2,
  FloatType const& twin_fraction)
{
  MMTBX_ASSERT(fo.size()==fc1.size());
  MMTBX_ASSERT(fo.size()==fc2.size());
  af::shared<FloatType> fc_abs(fo.size());
  for(std::size_t i=0; i < fo.size(); i++) {
    FloatType fc_abs1 = std::abs(fc1[i]);
    FloatType fc_abs2 = std::abs(fc2[i]);
    fc_abs[i]=std::sqrt((1-twin_fraction)*fc_abs1*fc_abs1+
                           twin_fraction *fc_abs2*fc_abs2);
  }
  return scale(fo,fc_abs.const_ref());
};

template <typename FloatType>
FloatType
r_factor(
  af::const_ref<FloatType> const& fo,
  af::const_ref<FloatType> const& fc,
  FloatType const& scale)
{
  MMTBX_ASSERT(fo.size()==fc.size());
  FloatType num=0.0;
  FloatType denum=0.0;
  for(std::size_t i=0; i < fo.size(); i++) {
    num += std::abs(fo[i] - fc[i] * scale);
    denum += fo[i];
  }
  if(denum == 0) return 1.e+9;
  return num/denum;
};

template <typename FloatType>
FloatType
r_factor(
  af::const_ref<FloatType> const& fo,
  af::const_ref< std::complex<FloatType> > const& fc,
  FloatType const& scale)
{
  MMTBX_ASSERT(fo.size()==fc.size());
  FloatType num=0.0;
  FloatType denum=0.0;
  for(std::size_t i=0; i < fo.size(); i++) {
    num += std::abs(fo[i] - std::abs(fc[i]) * scale);
    denum += fo[i];
  }
  if(denum == 0) return 1.e+9;
  return num/denum;
};

template <typename FloatType>
FloatType
r_factor(
  af::const_ref<FloatType> const& fo,
  af::const_ref<FloatType> const& fc)
{
  MMTBX_ASSERT(fo.size()==fc.size());
  FloatType sc = scale(fo,fc);
  return r_factor(fo,fc,sc);
};

template <typename FloatType, typename ComplexType>
FloatType
r_factor(
  af::const_ref<FloatType> const& fo,
  af::const_ref<std::complex<ComplexType> > const& fc)
{
  MMTBX_ASSERT(fo.size()==fc.size());
  FloatType sc = scale(fo,fc);
  return r_factor(fo,fc,sc);
};

template <typename FloatType, typename ComplexType>
FloatType
r_factor(
  af::const_ref<FloatType> const& fo,
  af::const_ref< std::complex<ComplexType> > const& fc1,
  af::const_ref< std::complex<ComplexType> > const& fc2,
  FloatType const& twin_fraction)
{
  MMTBX_ASSERT(fo.size()==fc1.size());
  MMTBX_ASSERT(fo.size()==fc2.size());
  af::shared<FloatType> fc_abs(fo.size());
  for(std::size_t i=0; i < fo.size(); i++) {
    FloatType fc_abs1 = std::abs(fc1[i]);
    FloatType fc_abs2 = std::abs(fc2[i]);
    fc_abs[i]=std::sqrt((1-twin_fraction)*fc_abs1*fc_abs1+
                           twin_fraction *fc_abs2*fc_abs2);
  }
  FloatType sc = scale(fo,fc_abs.const_ref());
  return r_factor(fo,fc_abs.const_ref(),sc);
};

template <typename FloatType, typename ComplexType>
FloatType
r_factor(
  af::const_ref<FloatType> const& fo,
  af::const_ref< std::complex<ComplexType> > const& fc1,
  af::const_ref< std::complex<ComplexType> > const& fc2,
  FloatType const& twin_fraction,
  FloatType const& scale)
{
  MMTBX_ASSERT(fo.size()==fc1.size());
  MMTBX_ASSERT(fo.size()==fc2.size());
  af::shared<FloatType> fc_abs(fo.size());
  for(std::size_t i=0; i < fo.size(); i++) {
    FloatType fc_abs1 = std::abs(fc1[i]);
    FloatType fc_abs2 = std::abs(fc2[i]);
    fc_abs[i]=std::sqrt((1-twin_fraction)*fc_abs1*fc_abs1+
                           twin_fraction *fc_abs2*fc_abs2);
  }
  return r_factor(fo,fc_abs.const_ref(),scale);
};
//------------------------------------------------------------------------------
template <typename FloatType, typename ComplexType>
vec3<FloatType>
ksol_bsol_grid_search(
  af::const_ref<FloatType> const& fo,
  af::const_ref<ComplexType> const& fc,
  af::const_ref<ComplexType> const& fm,
  sym_mat3<FloatType> const& b_cart,
  af::const_ref<FloatType> const& ksol_range,
  af::const_ref<FloatType> const& bsol_range,
  FloatType const& r_ref,
  af::const_ref<cctbx::miller::index<> > const& hkl,
  cctbx::uctbx::unit_cell const& uc)
{
  MMTBX_ASSERT(hkl.size() == fo.size());
  MMTBX_ASSERT(fo.size() == fc.size() && fc.size() == fm.size());
  af::shared<FloatType> s_mem = uc.stol_sq(hkl);
  af::const_ref<FloatType> ss = s_mem.const_ref();
  mat3<FloatType> a = uc.fractionalization_matrix();
  sym_mat3<FloatType> u_star = sym_mat3<FloatType> (b_cart).tensor_transform(a);
  FloatType k_best = 0.0;
  FloatType b_best = 0.0;
  FloatType r_best = r_ref;
  for(std::size_t i=0; i < ksol_range.size(); i++) {
    for(std::size_t j=0; j < bsol_range.size(); j++) {
      af::shared<ComplexType> f_model = f_model::f_model(
        fc, fm, u_star, ksol_range[i],bsol_range[j], ss, hkl);
      FloatType r = r_factor(fo, f_model.const_ref());
      if(r < r_best) {
        k_best = ksol_range[i];
        b_best = bsol_range[j];
        r_best = r;
      }
    }
  }
  return vec3<FloatType> (k_best,b_best,r_best);
};

template <typename FloatType, typename ComplexType>
vec3<FloatType>
ksol_bsol_grid_search(
  af::const_ref<FloatType> const& fo,
  af::const_ref<ComplexType> const& fc1,
  af::const_ref<ComplexType> const& fc2,
  af::const_ref<ComplexType> const& fm1,
  af::const_ref<ComplexType> const& fm2,
  sym_mat3<FloatType> const& b_cart,
  af::const_ref<FloatType> const& ksol_range,
  af::const_ref<FloatType> const& bsol_range,
  FloatType const& r_ref,
  af::const_ref<cctbx::miller::index<> > const& hkl,
  cctbx::uctbx::unit_cell const& uc,
  FloatType const& twin_fraction)
{
  MMTBX_ASSERT(hkl.size() == fo.size());
  MMTBX_ASSERT(fo.size() == fc1.size() && fc1.size() == fm1.size());
  MMTBX_ASSERT(fo.size() == fc2.size() && fc2.size() == fm2.size());
  af::shared<FloatType> s_mem = uc.stol_sq(hkl);
  af::const_ref<FloatType> ss = s_mem.const_ref();
  mat3<FloatType> a = uc.fractionalization_matrix();
  sym_mat3<FloatType> u_star = sym_mat3<FloatType> (b_cart).tensor_transform(a);
  FloatType k_best = 0.0;
  FloatType b_best = 0.0;
  FloatType r_best = r_ref;
  for(std::size_t i=0; i < ksol_range.size(); i++) {
    for(std::size_t j=0; j < bsol_range.size(); j++) {
      af::shared<ComplexType> f_model1 = f_model::f_model(
        fc1, fm1, u_star, ksol_range[i],bsol_range[j], ss, hkl);
      af::shared<ComplexType> f_model2 = f_model::f_model(
        fc2, fm2, u_star, ksol_range[i],bsol_range[j], ss, hkl);
      FloatType r = r_factor(fo, f_model1.const_ref(), f_model2.const_ref(),
        twin_fraction);
      if(r < r_best) {
        k_best = ksol_range[i];
        b_best = bsol_range[j];
        r_best = r;
      }
    }
  }
  return vec3<FloatType> (k_best,b_best,r_best);
};

////////////////////////////////////////////////////////////////////////////////
template <typename DataType, typename TagType>
void
  symmetrize_mask(
      af::ref<DataType, af::c_grid<3> > const& data,
      af::const_ref<TagType, af::c_grid<3> > const& tags)
  {
    MMTBX_ASSERT(tags.accessor().all_eq(data.accessor()));
    for(std::size_t i=0;i<data.size();i++) {
      if (tags[i] < 0) continue;
      if (data[i] == 0) data[tags[i]] = 0;
    }
    for(std::size_t i=0;i<data.size();i++) {
      if (tags[i] < 0) continue;
      data[i] = data[tags[i]];
    }
  }



af::shared<double> fb_cart(sym_mat3<double> const& b_cart,
                            af::const_ref<cctbx::miller::index<> > const& hkl,
                            cctbx::uctbx::unit_cell const& uc);



class target_gradients_aniso_ml {
public:
   target_gradients_aniso_ml(af::const_ref<double> const& fo,
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
                            double k);
   double target() const { return tgx; }
   af::shared<double> grad_b_cart() { return gtgx_u; }
   double grad_ksol() const { return gtgx_ksol; }
   double grad_bsol() const { return gtgx_bsol; }
   double grad_k() const { return gtgx_k; }
private:
   double tgx, gtgx_ksol, gtgx_bsol, gtgx_k;
   af::shared<double> gtgx_u;
};

}} // namespace mmtbx::bulk_solvent

#endif // MMTBX_BULK_SOLVENT_BULK_SOLVENT_H
