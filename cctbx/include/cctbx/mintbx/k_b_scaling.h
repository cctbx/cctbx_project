/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MINTBX_K_B_SCALING_H
#define CCTBX_MINTBX_K_B_SCALING_H

#include <cctbx/adptbx.h>
#include <scitbx/array_family/misc_functions.h>

namespace cctbx { namespace mintbx {

  template <typename FloatType = double>
  class k_b_scaling_target_and_gradients
  {
    public:
      typedef FloatType float_type;

      k_b_scaling_target_and_gradients() {}

      k_b_scaling_target_and_gradients(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<int> const& multiplicities,
        af::const_ref<FloatType> const& data_reference,
        af::const_ref<FloatType> const& data_scaled,
        FloatType k,
        FloatType b_iso,
        bool calc_gradient_k,
        bool calc_gradient_b)
      :
        anisotropic_flag_(false)
      {
        calculate(
          unit_cell, miller_indices, multiplicities,
          data_reference, data_scaled,
          k, b_iso, scitbx::sym_mat3<FloatType>(),
          calc_gradient_k, calc_gradient_b);
      }

      k_b_scaling_target_and_gradients(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<int> const& multiplicities,
        af::const_ref<FloatType> const& data_reference,
        af::const_ref<FloatType> const& data_scaled,
        FloatType k,
        scitbx::sym_mat3<FloatType> const& b_cif,
        bool calc_gradient_k,
        bool calc_gradient_b)
      :
        anisotropic_flag_(true)
      {
        calculate(
          unit_cell, miller_indices, multiplicities,
          data_reference, data_scaled,
          k, 0, b_cif,
          calc_gradient_k, calc_gradient_b);
      }

      FloatType
      target() const { return target_; }

      FloatType
      gradient_k() const { return gradient_k_; }

      bool
      anisotropic_flag() const { return anisotropic_flag_; }

      FloatType
      gradient_b_iso() const { return gradient_b_iso_; }

      scitbx::sym_mat3<FloatType> const&
      gradients_b_cif() const { return gradients_b_cif_; }

    protected:
      /* Mathematica:
           f=(d_sca - d_ref * k * Exp[c0*b0+c1*b1+c2*b2])^2
           D[f,k]
           D[f,b0]
       */
      void calculate(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<int> const& multiplicities,
        af::const_ref<FloatType> const& data_reference,
        af::const_ref<FloatType> const& data_scaled,
        FloatType k,
        FloatType b_iso,
        scitbx::sym_mat3<FloatType> const& b_cif,
        bool calc_gradient_k,
        bool calc_gradient_b);

      bool anisotropic_flag_;
      FloatType target_;
      FloatType gradient_k_;
      FloatType gradient_b_iso_;
      scitbx::sym_mat3<FloatType> gradients_b_cif_;
  };

  /* Mathematica input:
       f=(d_sca - d_ref * k * Exp[c0*b0+c1*b1+c2*b2])^2
       D[f,k]
       D[f,b0]
   */
  template <typename FloatType>
  void k_b_scaling_target_and_gradients<FloatType>::calculate(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<miller::index<> > const& miller_indices,
    af::const_ref<int> const& multiplicities,
    af::const_ref<FloatType> const& data_reference,
    af::const_ref<FloatType> const& data_scaled,
    FloatType k,
    FloatType b_iso,
    scitbx::sym_mat3<FloatType> const& b_cif,
    bool calc_gradient_k,
    bool calc_gradient_b)
  {
    if (multiplicities.size()) {
      CCTBX_ASSERT(miller_indices.size() == multiplicities.size());
    }
    CCTBX_ASSERT(miller_indices.size() == data_reference.size());
    CCTBX_ASSERT(miller_indices.size() == data_scaled.size());
    target_ = 0;
    gradient_k_ = 0;
    gradient_b_iso_ = 0;
    gradients_b_cif_.fill(0);
    FloatType sum_mult_data_scaled_sq = 0;
    FloatType mult = 1;
    std::size_t i;
    for(i=0;i<miller_indices.size();i++) {
      if (multiplicities.size()) mult = multiplicities[i];
      sum_mult_data_scaled_sq += mult * scitbx::fn::pow2(data_scaled[i]);
    }
    if (sum_mult_data_scaled_sq < 1.) sum_mult_data_scaled_sq = 1.;
    scitbx::sym_mat3<FloatType> u_star;
    if (anisotropic_flag_) {
      u_star = adptbx::u_cif_as_u_star(unit_cell, adptbx::b_as_u(b_cif));
    }
    for(i=0;i<miller_indices.size();i++) {
      if (multiplicities.size()) mult = multiplicities[i];
      FloatType dw, stol_sq;
      if (anisotropic_flag_) {
        dw = adptbx::debye_waller_factor_u_star(miller_indices[i], u_star);
      }
      else {
        stol_sq = unit_cell.stol_sq(miller_indices[i]);
        dw = adptbx::debye_waller_factor_b_iso(stol_sq, b_iso);
      }
      FloatType drkdw = data_reference[i] * k * dw;
      FloatType diff = data_scaled[i] - drkdw;
      target_ += mult * scitbx::fn::pow2(
        diff / std::sqrt(sum_mult_data_scaled_sq));
      if (calc_gradient_k || calc_gradient_b) {
        FloatType gk = -2 * mult * drkdw * diff / sum_mult_data_scaled_sq;
        if (calc_gradient_k) gradient_k_ += gk;
        if (calc_gradient_b) {
          if (anisotropic_flag_) {
            scitbx::sym_mat3<FloatType> dwc
              = adptbx::b_as_u(adptbx::u_cif_as_u_star(unit_cell,
                  adptbx::debye_waller_factor_u_star_coefficients(
                    miller_indices[i], scitbx::type_holder<FloatType>())));
            for(std::size_t j=0;j<dwc.size();j++) {
              gradients_b_cif_[j] += dwc[j] * gk;
            }
          }
          else {
            gradient_b_iso_ -= stol_sq * gk;
          }
        }
      }
    }
    if (calc_gradient_b && anisotropic_flag_) {
      gradients_b_cif_ *= -scitbx::constants::two_pi_sq;
    }
  }

}} // namespace cctbx::mintbx

#endif // CCTBX_MINTBX_K_B_SCALING_H
