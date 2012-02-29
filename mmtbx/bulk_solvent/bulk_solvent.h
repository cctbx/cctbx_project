#ifndef MMTBX_BULK_SOLVENT_BULK_SOLVENT_H
#define MMTBX_BULK_SOLVENT_BULK_SOLVENT_H

#include <mmtbx/error.h>
#include <mmtbx/import_scitbx_af.h>
#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <mmtbx/f_model/f_model.h>
#include <scitbx/matrix/outer_product.h>
#include <scitbx/array_family/versa_algebra.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/small_algebra.h>
#include <scitbx/matrix/eigensystem.h>
#include <scitbx/matrix/packed.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/math/cubic_equation.h>

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
// All scales (overall, overall anisotropic, etc) must be applied.
template <
  typename FloatType=double,
  typename ComplexType=std::complex<double>,
  typename CubicEqType=scitbx::math::cubic_equation::real<long double, double> >
class bulk_solvent_scale_coefficients_analytical
{
public:
  af::shared<FloatType> x, r;
  FloatType x_best, r_best;

  bulk_solvent_scale_coefficients_analytical() {}

  bulk_solvent_scale_coefficients_analytical(
    af::const_ref<FloatType> const& f_obs,
    af::const_ref<ComplexType> const& f_calc,
    af::const_ref<ComplexType> const& f_mask,
    af::const_ref<bool> const& selection)
  :
  x_best(0), r_best(0)
  {
    MMTBX_ASSERT(f_obs.size() == f_calc.size());
    MMTBX_ASSERT(f_obs.size() == f_mask.size());
    MMTBX_ASSERT(f_obs.size() == selection.size());
    FloatType a2 = 0.0;
    FloatType d3 = 0.0;
    FloatType a = 0.0;
    FloatType b = 0.0;
    FloatType c = 0.0;
    for(std::size_t i=0; i < f_obs.size(); i++) {
      if(selection[i]) {
        FloatType p = std::real(f_calc[i]);
        FloatType r = std::imag(f_calc[i]);
        FloatType q = std::real(f_mask[i]);
        FloatType t = std::imag(f_mask[i]);
        FloatType I = f_obs[i]*f_obs[i];
        FloatType v = p*q+r*t;
        FloatType w = q*q+t*t;
        FloatType u = p*p+r*r;
        a2 += (u*I);
        d3 += (w*w);
        a  += (3.*v*w);
        b  += (2*v*v + u*w - w*I);
        c  += (u*v - v*I);
      }
    }
    MMTBX_ASSERT(d3 != 0.0);
    // coefficients of x**3 + ax**2 + bx + c = 0
    CubicEqType ceo = CubicEqType(1,a/d3,b/d3,c/d3);
    // we are interested in non-negative roots only
    x.push_back(0);
    for(std::size_t j=0; j < 3; j++) {
      if(ceo.x[j]) {
        FloatType root = *ceo.x[j];
        if(root>=0) {
          x.push_back(root);
        }
        MMTBX_ASSERT(std::abs(*ceo.residual()[j]) < 1.e-4);
      }
    }
    // put together plausible results
    bool zero = false;
    af::shared<ComplexType> f_model(f_obs.size());
    for(std::size_t j=0; j < x.size(); j++) {
      if((x[j]>0) || (x[j]==0 && zero==false)) {
        for(std::size_t i=0; i < f_obs.size(); i++) {
          if(selection[i]) {
            f_model[i] = f_calc[i] + x[j] * f_mask[i];
          }
        }
        r.push_back(r_factor(f_obs, f_model.const_ref()));
      }
      else {
        r.push_back(-1);
      }
      if(x[j]==0) zero = true;
    }
    // select best result
    x_best = x[0];
    r_best=1.e+9;
    for(std::size_t j=0; j < x.size(); j++) {
      if(r[j]<=r_best && r[j]>=0) {
        r_best = r[j];
        x_best = x[j];
      }
    }
  }
};
//------------------------------------------------------------------------------

template <
  typename FloatType=double,
  typename ComplexType=std::complex<double>,
  typename CubicEqType=scitbx::math::cubic_equation::real<long double,double> >
class overall_and_bulk_solvent_scale_coefficients_analytical
{
public:
  af::shared<FloatType> x, r;
  FloatType x_best, r_best;

  overall_and_bulk_solvent_scale_coefficients_analytical() {}

  overall_and_bulk_solvent_scale_coefficients_analytical(
    af::const_ref<FloatType> const& f_obs,
    af::const_ref<ComplexType> const& f_calc,
    af::const_ref<ComplexType> const& f_mask,
    af::const_ref<bool> const& selection)
  :
  x_best(0), r_best(0)
  {
    MMTBX_ASSERT(f_obs.size() == f_calc.size());
    MMTBX_ASSERT(f_obs.size() == f_mask.size());
    MMTBX_ASSERT(f_obs.size() == selection.size());
    FloatType a2 = 0.0;
    FloatType b2 = 0.0;
    FloatType c2 = 0.0;
    FloatType y2 = 0.0;
    FloatType a3 = 0.0;
    FloatType b3 = 0.0;
    FloatType c3 = 0.0;
    FloatType d3 = 0.0;
    FloatType y3 = 0.0;
    for(std::size_t i=0; i < f_obs.size(); i++) {
      if(selection[i]) {
        FloatType p = std::real(f_calc[i]);
        FloatType r = std::imag(f_calc[i]);
        FloatType q = std::real(f_mask[i]);
        FloatType t = std::imag(f_mask[i]);
        FloatType I = f_obs[i]*f_obs[i];
        FloatType v = p*q+r*t;
        FloatType w = q*q+t*t;
        FloatType u = p*p+r*r;
        a2 += (u*I);
        b2 += (2.*v*I);
        c2 += (w*I);
        y2 += (I*I);
        a3 += (u*v);
        b3 += (2.*v*v+u*w);
        c3 += (3.*v*w);
        d3 += (w*w);
        y3 += (v*I);
      }
    }
    FloatType den = d3*y2-c2*c2;
    MMTBX_ASSERT(den != 0.0);
    // coefficients of x**3 + ax**2 + bc + c = 0
    FloatType a = (c3*y2-c2*b2-c2*y3)/den;
    FloatType b = (b3*y2-c2*a2-y3*b2)/den;
    FloatType c = (a3*y2-y3*a2)/den;
    CubicEqType ceo = CubicEqType(1,a,b,c);
    // we are interested in non-negative roots only
    x.push_back(0);
    for(std::size_t j=0; j < 3; j++) {
      if(ceo.x[j]) {
        FloatType root = *ceo.x[j];
        if(root>=0) {
          x.push_back(root);
        }
        MMTBX_ASSERT(std::abs(*ceo.residual()[j]) < 1.e-4);
      }
    }
    // put together plausible results
    bool zero = false;
    af::shared<ComplexType> f_model(f_obs.size());
    for(std::size_t j=0; j < x.size(); j++) {
      if((x[j]>0) || (x[j]==0 && zero==false)) {
        for(std::size_t i=0; i < f_obs.size(); i++) {
          if(selection[i]) {
            f_model[i] = f_calc[i] + x[j] * f_mask[i];
          }
        }
        r.push_back(r_factor(f_obs, f_model.const_ref()));
      }
      else {
        r.push_back(-1);
      }
      if(x[j]==0) zero = true;
    }
    // select best result
    x_best = x[0];
    r_best=1.e+9;
    for(std::size_t j=0; j < x.size(); j++) {
      if(r[j]<=r_best && r[j]>=0) {
        r_best = r[j];
        x_best = x[j];
      }
    }
  }
};

template <typename FloatType=double,
          typename ComplexType=std::complex<double> >
class aniso_u_scaler
{
public:
  std::size_t n_rows;
  af::shared<FloatType> u_star_independent;
  af::shared<FloatType> a;

  aniso_u_scaler() {}

  aniso_u_scaler(
    af::const_ref<ComplexType> const& f_model,
    af::const_ref<FloatType> const& f_obs,
    af::const_ref<cctbx::miller::index<> > const& miller_indices,
    af::const_ref<FloatType, af::mat_grid> const& adp_constraint_matrix)
  :
  n_rows(adp_constraint_matrix.accessor().n_rows()),
  u_star_independent(n_rows, 0)
  {
    MMTBX_ASSERT(f_obs.size() == f_model.size());
    MMTBX_ASSERT(f_obs.size() == miller_indices.size());
    FloatType minus_two_pi_sq = -2.*std::pow(scitbx::constants::pi, 2);
    af::versa<FloatType, af::mat_grid> m_(af::mat_grid(n_rows, n_rows), 0);
    af::versa<FloatType, af::mat_grid> m(af::mat_grid(n_rows, n_rows), 0);
    af::small<FloatType, 6> b(n_rows, 0);
    af::small<FloatType, 6> vr(n_rows);
    for(std::size_t i=0; i < f_obs.size(); i++) {
      cctbx::miller::index<> const& miller_index = miller_indices[i];
      int i0=miller_index[0],i1=miller_index[1],i2=miller_index[2];
      FloatType fm_abs = std::abs(f_model[i]);
      FloatType fo_i = f_obs[i];
      MMTBX_ASSERT(fm_abs > 0);
      MMTBX_ASSERT(fo_i > 0);
      FloatType z = std::log(fo_i/fm_abs)/minus_two_pi_sq;
#define _ static_cast<FloatType>
      FloatType const v[] = {
        _(i0*i0), _(i1*i1), _(i2*i2), _(2*i0*i1), _(2*i0*i2), _(2*i1*i2)};
#undef _
      scitbx::matrix::multiply(
        /*a*/ adp_constraint_matrix.begin(),
        /*b*/ v,
        /*ar*/ n_rows,
        /*ac*/ 6,
        /*bc*/ 1,
        /*ab*/ vr.begin());
      scitbx::matrix::outer_product(m_.begin(),vr.const_ref(),vr.const_ref());
      m += m_;
      b += z*vr;
    }
    af::versa<FloatType, af::c_grid<2> > m_inv(
       scitbx::matrix::packed_u_as_symmetric(
         scitbx::matrix::eigensystem::real_symmetric<FloatType>(
           m.const_ref(), /*relative_epsilon*/ 1.e-9,/*absolute_epsilon*/ 1.e-9)
             .generalized_inverse_as_packed_u().const_ref()));
    u_star_independent = af::matrix_multiply(m_inv.const_ref(), b.const_ref());
  }

  /* SHELXL-like scale: Acta Cryst. (1999). D55, 1158-1167, page 1163
     Note: in the paper it is applied to Fcalc^2, here it is applied to Fcalc.
           The code computes coefficients a. */
  /* XXX try 'long double' to see if that pushed regression test tolerances
     to lower values */
  aniso_u_scaler(
    af::const_ref<ComplexType> const& f_model,
    af::const_ref<FloatType> const& f_obs,
    af::const_ref<cctbx::miller::index<> > const& miller_indices,
    cctbx::uctbx::unit_cell const& unit_cell)
  :
  a(12, 0)
  {
    MMTBX_ASSERT(f_obs.size() == f_model.size());
    MMTBX_ASSERT(f_obs.size() == miller_indices.size());
    af::versa<FloatType, af::mat_grid> m_(af::mat_grid(12, 12), 0);
    af::versa<FloatType, af::mat_grid> m(af::mat_grid(12, 12), 0);
    af::tiny<FloatType, 12> b;
    b.fill(0);
    af::double6 p = unit_cell.reciprocal_parameters();
    FloatType as=p[0], bs=p[1], cs=p[2];
    for(std::size_t i=0; i < f_obs.size(); i++) {
      cctbx::miller::index<> const& miller_index = miller_indices[i];
      int h=miller_index[0], k=miller_index[1], l=miller_index[2];
      FloatType fm_i = std::abs(f_model[i]);
      FloatType stol = unit_cell.stol_sq(miller_index);
      FloatType s = 0;
      if(stol != 0) s = 1./stol;
      af::tiny<FloatType, 12> v_;
      FloatType* v = v_.begin();
      *v++ = h*h*as*as*s;
      *v++ = h*h*as*as;
      *v++ = k*k*bs*bs*s;
      *v++ = k*k*bs*bs;
      *v++ = l*l*cs*cs*s;
      *v++ = l*l*cs*cs;
      *v++ = 2*k*l*bs*cs*s;
      *v++ = 2*k*l*bs*cs;
      *v++ = 2*h*l*as*cs*s;
      *v++ = 2*h*l*as*cs;
      *v++ = 2*h*k*as*bs*s;
      *v++ = 2*h*k*as*bs;
      b += (f_obs[i]-fm_i)*fm_i*v_;
      v_ *= fm_i;
      scitbx::matrix::outer_product(m_.begin(),v_.const_ref(),v_.const_ref());
      m += m_;
    }
    af::versa<FloatType, af::c_grid<2> > m_inv(
       scitbx::matrix::packed_u_as_symmetric(
         scitbx::matrix::eigensystem::real_symmetric<FloatType>(
           m.const_ref(), /*relative_epsilon*/ 1.e-9,/*absolute_epsilon*/ 1.e-9)
             .generalized_inverse_as_packed_u().const_ref()));
    a = af::matrix_multiply(m_inv.const_ref(), b.const_ref());
  }
};

template <typename FloatType, typename ComplexType>
 af::tiny<FloatType, 2>
 k_mask_and_k_overall_grid_search(
   af::const_ref<FloatType>   const& f_obs,
   af::const_ref<ComplexType> const& f_calc,
   af::const_ref<ComplexType> const& f_mask,
   af::const_ref<FloatType>   const& k_mask_range
   )
 {
   MMTBX_ASSERT(f_mask.size() == f_obs.size());
   MMTBX_ASSERT(f_obs.size() == f_calc.size());
   FloatType k_mask_best = 0.0;
   FloatType k_overall_best = 0.0;
   FloatType r_best = r_factor(f_obs, f_calc);
   af::shared<ComplexType> f_model(f_obs.size());
   for(std::size_t i=0; i < k_mask_range.size(); i++) {
     FloatType k_mask = k_mask_range[i];
     for(std::size_t j=0; j < f_obs.size(); j++) {
       f_model[j] = f_calc[j] + k_mask * f_mask[j];
     }
     FloatType k_overall = scale(f_obs, f_model.const_ref());
     FloatType r = r_factor(f_obs, f_model.const_ref(), k_overall);
     if(r < r_best) {
       k_mask_best = k_mask;
       k_overall_best = k_overall;
       r_best = r;
     }
   }
   return af::tiny<FloatType, 2> (k_mask_best, k_overall_best);
 };

template <typename FloatType, typename ComplexType>
 af::shared<FloatType>
 ksol_bsol_grid_search(
   af::const_ref<FloatType>   const& f_obs,
   af::const_ref<ComplexType> const& f_calc,
   af::const_ref<ComplexType> const& f_mask,
   af::const_ref<FloatType>   const& k_sol_range,
   af::const_ref<FloatType>   const& b_sol_range,
   af::const_ref<FloatType>   const& ss,
   FloatType                  const& scalar_scale,
   af::const_ref<FloatType>   const& overall_scale,
   af::const_ref<FloatType>   const& overall_anisotropic_scale,
   FloatType                  const& r_ref)
 {
   MMTBX_ASSERT(f_mask.size() == f_obs.size());
   MMTBX_ASSERT(f_obs.size() == f_calc.size());
   MMTBX_ASSERT(ss.size() == f_calc.size());
   MMTBX_ASSERT(overall_scale.size() == f_calc.size());
   MMTBX_ASSERT(overall_anisotropic_scale.size() == f_calc.size());
   FloatType k_best = 0.0;
   FloatType b_best = 0.0;
   FloatType r_best = r_ref;
   af::shared<ComplexType> f_model(ss.size());
   af::shared<FloatType> bulk_solvent_scale(f_obs.size());
   for(std::size_t i=0; i < k_sol_range.size(); i++) {
     FloatType ks = k_sol_range[i];
     for(std::size_t j=0; j < b_sol_range.size(); j++) {
       FloatType mbs = -b_sol_range[j];
       for(std::size_t k=0; k < f_obs.size(); k++) {
         FloatType kbs = ks * std::exp(mbs * ss[k]);
         f_model[k] = scalar_scale * overall_scale[k] *
           overall_anisotropic_scale[k] * (f_calc[k] + kbs * f_mask[k]);
       }
       FloatType r = r_factor(f_obs, f_model.const_ref());
       if(r < r_best) {
         k_best = k_sol_range[i];
         b_best = b_sol_range[j];
         r_best = r;
       }
     }
   }
   for(std::size_t k=0; k < f_obs.size(); k++) {
     bulk_solvent_scale[k] = k_best * std::exp(-b_best * ss[k]);
   }
   return bulk_solvent_scale;
 };

template <typename FloatType>
 af::shared<FloatType>
 set_to_liear_interpolated(
   af::const_ref<FloatType> const& ss,
   FloatType                const& k,
   FloatType                const& b,
   af::const_ref<bool>      const& selection,
   af::shared<FloatType>           data)
 {
   af::shared<FloatType> k_mask(ss.size());
   for(std::size_t i=0; i < ss.size(); i++) {
     if(selection[i]) {
       FloatType v = k * ss[i] + b;
       if(v<0) v=0;
       data[i] = v;
     }
   }
   return data;
 };

template <typename FloatType>
 af::shared<FloatType>
 set_k_mask_to_cubic_polynom(
   af::const_ref<FloatType> const& ss,
   FloatType                const& ss_cutoff,
   af::tiny<FloatType, 4>   const& coeffs)
 {
   af::shared<FloatType> k_mask(ss.size());
   for(std::size_t i=0; i < ss.size(); i++) {
     FloatType x = ss[i];
     FloatType f = coeffs[0] + coeffs[1]*x + coeffs[2]*x*x + coeffs[3]*x*x*x;
     if(f<0) f = 0;
     if(x<ss_cutoff) k_mask[i] = f;
     else k_mask[i] = 0;
   }
   return k_mask;
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

template <typename FloatType, typename ComplexType>
FloatType
scale(af::const_ref<FloatType> const& fo,
      af::const_ref<ComplexType> const& fc,
      af::const_ref<bool> const& selection)
{
    MMTBX_ASSERT(fo.size()==fc.size());
    MMTBX_ASSERT(fo.size()==selection.size());
    FloatType num=0.0;
    FloatType denum=0.0;
    for(std::size_t i=0; i < fo.size(); i++) {
      if(selection[i]) {
        FloatType fc_abs = std::abs(fc[i]);
        num += fo[i] * fc_abs;
        denum += fc_abs * fc_abs;
      }
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
  af::const_ref< std::complex<FloatType> > const& fc,
  af::const_ref<bool> const& selection,
  FloatType const& scale)
{
  MMTBX_ASSERT(fo.size()==fc.size());
  MMTBX_ASSERT(fo.size()==selection.size());
  FloatType num=0.0;
  FloatType denum=0.0;
  for(std::size_t i=0; i < fo.size(); i++) {
    if(selection[i]) {
      num += std::abs(fo[i] - std::abs(fc[i]) * scale);
      denum += fo[i];
    }
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
  af::const_ref<std::complex<ComplexType> > const& fc,
  af::const_ref<bool> const& selection)
{
  MMTBX_ASSERT(fo.size()==fc.size());
  MMTBX_ASSERT(fo.size()==selection.size());
  FloatType sc = scale(fo,fc,selection);
  return r_factor(fo,fc,selection,sc);
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
