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

  class k_b_scaling_target_and_gradients
  {
    public:
      k_b_scaling_target_and_gradients() {}

      /* Anisotropic k scaling and B scaling

         Mathematica:

         f=(d_sca - d_ref * k * Exp[c0*b0+c1*b1+c2*b2])^2
         D[f,b0]

                   c0 b0 + c1 b1 + c2 b2
            -2 c0 E                      k d_ref

             c0 b0 + c1 b1 + c2 b2
         (-(E                      k d_ref) + d_sca)
       */
      k_b_scaling_target_and_gradients(
        uctbx::UnitCell const& ucell,
        af::shared<Miller::Index> miller_indices,
        af::shared<int> multiplicities,
        af::shared<double> data_reference,
        af::shared<double> data_scaled,
        double k,
        af::tiny<double, 6> const& b_cif,
        bool calc_gradient_k,
        bool calc_gradients_b_cif)
      {
        using namespace adptbx;
        if (multiplicities.size()) {
          cctbx_assert(miller_indices.size() == multiplicities.size());
        }
        cctbx_assert(miller_indices.size() == data_reference.size());
        cctbx_assert(miller_indices.size() == data_scaled.size());
        target_ = 0;
        if (calc_gradient_k) gradient_k_ = 0;
        if (calc_gradients_b_cif) gradients_b_cif_.fill(0);
        double sum_mult_data_scaled_sq = 0;
        double mult = 1;
        for(std::size_t i=0;i<miller_indices.size();i++) {
          if (multiplicities.size()) mult = multiplicities[i];
          sum_mult_data_scaled_sq += mult * math::pow2(data_scaled[i]);
        }
        if (sum_mult_data_scaled_sq < 1.) sum_mult_data_scaled_sq = 1.;
        af::tiny<double, 6> u_star = Ucif_as_Ustar(ucell, B_as_U(b_cif));
        for(std::size_t i=0;i<miller_indices.size();i++) {
          if (multiplicities.size()) mult = multiplicities[i];
          double dw = DebyeWallerFactorUstar(miller_indices[i], u_star);
          double drkdw = data_reference[i] * k * dw;
          double diff = data_scaled[i] - drkdw;
          target_ += mult * math::pow2(
            diff / std::sqrt(sum_mult_data_scaled_sq));
          if (calc_gradient_k || calc_gradients_b_cif) {
            double gk = -2 * mult * drkdw * diff / sum_mult_data_scaled_sq;
            if (calc_gradient_k) gradient_k_ += gk;
            if (calc_gradients_b_cif) {
              af::tiny<double, 6> dwc = B_as_U(Ucif_as_Ustar(ucell,
                DebyeWallerFactorUstarCoefficients(
                  miller_indices[i], type_holder<double>())));
              for(std::size_t j=0;j<dwc.size();j++) {
                gradients_b_cif_[j] += dwc[j] * gk;
              }
            }
          }
        }
      }

      double target() const { return target_; }

      double gradient_k() const { return gradient_k_; }

      af::tiny<double, 6> const& gradients_b_cif() const
      {
        return gradients_b_cif_;
      }

    protected:
      double target_;
      double gradient_k_;
      af::tiny<double, 6> gradients_b_cif_;
  };

}} // namespace cctbx::mintbx

#endif // CCTBX_MINTBX_K_B_SCALING_H
