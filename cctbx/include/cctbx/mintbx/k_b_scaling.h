// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MINTBX_K_B_SCALING_H
#define CCTBX_MINTBX_K_B_SCALING_H

#include <cctbx/adptbx.h>
#include <cctbx/sgtbx/groups.h>

namespace cctbx { namespace mintbx {

  template <typename FloatType>
  class k_b_scaling_target_and_gradients
  {
    public:
      k_b_scaling_target_and_gradients() {}

      k_b_scaling_target_and_gradients(
        uctbx::UnitCell const& ucell,
        af::shared<Miller::Index> miller_indices,
        af::shared<int> multiplicities,
        af::shared<FloatType> data_reference,
        af::shared<FloatType> data_scaled,
        FloatType k,
        FloatType b_iso,
        bool calc_gradient_k,
        bool calc_gradient_b)
      : anisotropic_(false)
      {
        calculate(
          ucell, miller_indices, multiplicities,
          data_reference, data_scaled,
          k, b_iso, af::tiny<FloatType, 6>(),
          calc_gradient_k, calc_gradient_b);
      }

      k_b_scaling_target_and_gradients(
        uctbx::UnitCell const& ucell,
        af::shared<Miller::Index> miller_indices,
        af::shared<int> multiplicities,
        af::shared<FloatType> data_reference,
        af::shared<FloatType> data_scaled,
        FloatType k,
        af::tiny<FloatType, 6> const& b_cif,
        bool calc_gradient_k,
        bool calc_gradient_b)
      : anisotropic_(true)
      {
        calculate(
          ucell, miller_indices, multiplicities,
          data_reference, data_scaled,
          k, 0, b_cif,
          calc_gradient_k, calc_gradient_b);
      }

      FloatType target() const { return target_; }

      FloatType gradient_k() const { return gradient_k_; }

      bool anisotropic() const { return anisotropic_; }

      FloatType gradient_b_iso() const { return gradient_b_iso_; }

      af::tiny<FloatType, 6> const& gradients_b_cif() const
      {
        return gradients_b_cif_;
      }

    protected:
      /* Mathematica:
           f=(d_sca - d_ref * k * Exp[c0*b0+c1*b1+c2*b2])^2
           D[f,k]
           D[f,b0]
       */
      void calculate(
        uctbx::UnitCell const& ucell,
        af::shared<Miller::Index> miller_indices,
        af::shared<int> multiplicities,
        af::shared<FloatType> data_reference,
        af::shared<FloatType> data_scaled,
        FloatType k,
        FloatType b_iso,
        af::tiny<FloatType, 6> const& b_cif,
        bool calc_gradient_k,
        bool calc_gradient_b);

      bool anisotropic_;
      FloatType target_;
      FloatType gradient_k_;
      FloatType gradient_b_iso_;
      af::tiny<FloatType, 6> gradients_b_cif_;
  };

  /* Mathematica input:
       f=(d_sca - d_ref * k * Exp[c0*b0+c1*b1+c2*b2])^2
       D[f,k]
       D[f,b0]
   */
  template <typename FloatType>
  void k_b_scaling_target_and_gradients<FloatType>::calculate(
    uctbx::UnitCell const& ucell,
    af::shared<Miller::Index> miller_indices,
    af::shared<int> multiplicities,
    af::shared<FloatType> data_reference,
    af::shared<FloatType> data_scaled,
    FloatType k,
    FloatType b_iso,
    af::tiny<FloatType, 6> const& b_cif,
    bool calc_gradient_k,
    bool calc_gradient_b)
  {
    using namespace adptbx;
    if (multiplicities.size()) {
      cctbx_assert(miller_indices.size() == multiplicities.size());
    }
    cctbx_assert(miller_indices.size() == data_reference.size());
    cctbx_assert(miller_indices.size() == data_scaled.size());
    target_ = 0;
    gradient_k_ = 0;
    gradient_b_iso_ = 0;
    gradients_b_cif_.fill(0);
    FloatType sum_mult_data_scaled_sq = 0;
    FloatType mult = 1;
    for(std::size_t i=0;i<miller_indices.size();i++) {
      if (multiplicities.size()) mult = multiplicities[i];
      sum_mult_data_scaled_sq += mult * math::pow2(data_scaled[i]);
    }
    if (sum_mult_data_scaled_sq < 1.) sum_mult_data_scaled_sq = 1.;
    af::tiny<FloatType, 6> u_star;
    if (anisotropic_) {
      u_star = Ucif_as_Ustar(ucell, B_as_U(b_cif));
    }
    for(std::size_t i=0;i<miller_indices.size();i++) {
      if (multiplicities.size()) mult = multiplicities[i];
      FloatType dw, stol2;
      if (anisotropic_) {
        dw = DebyeWallerFactorUstar(miller_indices[i], u_star);
      }
      else {
        stol2 = ucell.stol2(miller_indices[i]);
        dw = DebyeWallerFactorBiso(stol2, b_iso);
      }
      FloatType drkdw = data_reference[i] * k * dw;
      FloatType diff = data_scaled[i] - drkdw;
      target_ += mult * math::pow2(
        diff / std::sqrt(sum_mult_data_scaled_sq));
      if (calc_gradient_k || calc_gradient_b) {
        FloatType gk = -2 * mult * drkdw * diff / sum_mult_data_scaled_sq;
        if (calc_gradient_k) gradient_k_ += gk;
        if (calc_gradient_b) {
          if (anisotropic_) {
            af::tiny<FloatType, 6> dwc = B_as_U(Ucif_as_Ustar(ucell,
              DebyeWallerFactorUstarCoefficients(
                miller_indices[i], type_holder<FloatType>())));
            for(std::size_t j=0;j<dwc.size();j++) {
              gradients_b_cif_[j] += dwc[j] * gk;
            }
          }
          else {
            gradient_b_iso_ -= stol2 * gk;
          }
        }
      }
    }
  }

}} // namespace cctbx::mintbx

#endif // CCTBX_MINTBX_K_B_SCALING_H
