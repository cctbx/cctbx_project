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

         f=(d_sca - d_ref * k * Exp[c0*u0+c1*u1+c2*u2])^2
         D[f,u0]

                   c0 u0 + c1 u1 + c2 u2
            -2 c0 E                      k d_ref

             c0 u0 + c1 u1 + c2 u2
         (-(E                      k d_ref) + d_sca)
       */
      k_b_scaling_target_and_gradients(
        af::shared<Miller::Index> miller_indices,
        af::shared<int> multiplicities,
        af::shared<double> data_reference,
        af::shared<double> data_scaled,
        double k,
        af::tiny<double, 6> const& u_star,
        bool compute_gradients)
      {
        cctbx_assert(miller_indices.size() == multiplicities.size());
        cctbx_assert(miller_indices.size() == data_reference.size());
        cctbx_assert(miller_indices.size() == data_scaled.size());
        target_ = 0;
        if (compute_gradients) {
          gradient_k_ = 0;
          gradients_u_star_.fill(0);
        }
        double mult = 1;
        for(std::size_t i=0;i<miller_indices.size();i++) {
          if (multiplicities.size()) mult = multiplicities[i];
          double dw = adptbx::DebyeWallerFactorUstar(miller_indices[i], u_star);
          double drkdw  = mult * data_reference[i] * k * dw;
          double diff  = mult * data_scaled[i] - drkdw;
          target_ += math::pow2(diff);
          if (compute_gradients) {
            gradient_k_ = -2 * drkdw * diff;
            af::tiny<double, 6>
            dwc = adptbx::DebyeWallerFactorUstarCoefficients(miller_indices[i],
              type_holder<double>());
            for(std::size_t j=0;j<dwc.size();j++) {
              gradients_u_star_[j] += dwc[j] * gradient_k_;
            }
          }
        }
      }

      double target() const { return target_; }

      double gradient_k() const { return gradient_k_; }

      af::tiny<double, 6> const& gradients_u_star() const
      {
        return gradients_u_star_;
      }

    protected:
      double target_;
      double gradient_k_;
      af::tiny<double, 6> gradients_u_star_;
  };

}} // namespace cctbx::mintbx

#endif // CCTBX_MINTBX_K_B_SCALING_H
