// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Oct 2002: Modified copy of phenix/xray_targets.h (rwgk)
     Mar 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_XRAY_TARGETS_H
#define CCTBX_XRAY_TARGETS_H

#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <complex>
#include <cmath>

namespace cctbx { namespace xray { namespace targets {

  namespace detail {

    template <typename FobsValueType,
              typename WeightValueType,
              typename FcalcValueType>
    FobsValueType
    scale_factor_calculation(
      af::const_ref<FobsValueType> const& fobs,
      af::const_ref<WeightValueType> const& weights,
      af::const_ref<FcalcValueType> const& fcalc)
    {
      CCTBX_ASSERT(fobs.size() == weights.size() || weights.size() == 0);
      CCTBX_ASSERT(fobs.size() == fcalc.size());
      FobsValueType sum_w_fobs_fcalc(0);
      FobsValueType sum_w_fcalc2(0);
      if (weights.size()) {
        for(std::size_t i=0;i<fobs.size();i++) {
          FobsValueType abs_fcalc = std::abs(fcalc[i]);
          sum_w_fobs_fcalc += weights[i] * fobs[i] * abs_fcalc;
          sum_w_fcalc2 += weights[i] * abs_fcalc * abs_fcalc;
        }
      }
      else {
        for(std::size_t i=0;i<fobs.size();i++) {
          FobsValueType abs_fcalc = std::abs(fcalc[i]);
          sum_w_fobs_fcalc += fobs[i] * abs_fcalc;
          sum_w_fcalc2 += abs_fcalc * abs_fcalc;
        }
      }
      if (sum_w_fcalc2 == 0) {
        throw cctbx::error(
          "Cannot calculate scale factor: sum of weights * fcalc^2 == 0.");
      }
      return sum_w_fobs_fcalc / sum_w_fcalc2;
    }

    template <typename ValueValueType,
              typename WeightValueType>
    ValueValueType
    sum_weighted_values_squared(
      af::const_ref<ValueValueType> const& values,
      af::const_ref<WeightValueType> const& weights)
    {
      CCTBX_ASSERT(values.size() == weights.size() || weights.size() == 0);
      ValueValueType result = 0.;
      if (weights.size()) {
        for(std::size_t i=0;i<values.size();i++) {
          result += weights[i] * values[i] * values[i];
        }
      }
      else {
        for(std::size_t i=0;i<values.size();i++) {
          result += values[i] * values[i];
        }
      }
      return result;
    }

  } // namespace detail

  template <typename FobsValueType = double,
            typename WeightValueType = FobsValueType,
            typename FcalcValueType = std::complex<FobsValueType> >
  class least_squares_residual
  {
    public:
      least_squares_residual() {}

      least_squares_residual(
        af::const_ref<FobsValueType> const& fobs,
        af::const_ref<WeightValueType> const& weights,
        af::const_ref<FcalcValueType> const& fcalc,
        bool compute_derivatives=false,
        FobsValueType const& scale_factor=0)
      :
        scale_factor_(scale_factor)
      {
        init(fobs, weights, fcalc, compute_derivatives);
      }

      least_squares_residual(
        af::const_ref<FobsValueType> const& fobs,
        af::const_ref<FcalcValueType> const& fcalc,
        bool compute_derivatives=false,
        FobsValueType const& scale_factor=0)
      :
        scale_factor_(scale_factor)
      {
        init(fobs, af::const_ref<WeightValueType>(0,0),
             fcalc, compute_derivatives);
      }

      FobsValueType
      scale_factor() const { return scale_factor_; }

      FobsValueType
      target() const { return target_; }

      af::shared<FcalcValueType>
      derivatives() const
      {
        return derivatives_;
      }

    protected:
      FobsValueType scale_factor_;
      FobsValueType target_;
      af::shared<FcalcValueType> derivatives_;

      void init(
        af::const_ref<FobsValueType> const& fobs,
        af::const_ref<WeightValueType> const& weights,
        af::const_ref<FcalcValueType> const& fcalc,
        bool compute_derivatives);
  };

  template <typename FobsValueType,
            typename WeightValueType,
            typename FcalcValueType>
  void
  least_squares_residual<FobsValueType, WeightValueType, FcalcValueType>
  ::init(
    af::const_ref<FobsValueType> const& fobs,
    af::const_ref<WeightValueType> const& weights,
    af::const_ref<FcalcValueType> const& fcalc,
    bool compute_derivatives)
  {
    if (scale_factor_ == 0) {
      scale_factor_ = detail::scale_factor_calculation(
        fobs, weights, fcalc);
    }
    FobsValueType sum_w_fobs2 = detail::sum_weighted_values_squared(
      fobs, weights);
    if (sum_w_fobs2 == 0) {
      throw cctbx::error(
        "Cannot calculate least-squares residual:"
        " sum of weights * fobs^2 == 0.");
    }
    target_ = 0;
    if (compute_derivatives) {
      derivatives_ = af::shared<FcalcValueType>(fobs.size());
    }
    WeightValueType w(1);
    for(std::size_t i=0;i<fobs.size();i++) {
      FobsValueType abs_fcalc = std::abs(fcalc[i]);
      FobsValueType delta = fobs[i] - scale_factor_ * abs_fcalc;
      if (weights.size()) w = weights[i];
      target_ += w * delta * delta;
      if (compute_derivatives && abs_fcalc != 0) {
        derivatives_[i] = -2. * scale_factor_ * w * delta
                        / (sum_w_fobs2 * abs_fcalc) * std::conj(fcalc[i]);
      }
    }
    target_ /= sum_w_fobs2;
  }

  template <typename FobsValueType = double,
            typename WeightValueType = int,
            typename FcalcValueType = std::complex<FobsValueType>,
            typename SumWeightsType = long>
  class intensity_correlation
  {
    public:
      intensity_correlation() {}

      intensity_correlation(
        af::const_ref<FobsValueType> const& fobs,
        af::const_ref<WeightValueType> const& weights,
        af::const_ref<FcalcValueType> const& fcalc,
        bool compute_derivatives = false)
      {
        init(fobs, weights, fcalc, compute_derivatives);
      }

      intensity_correlation(
        af::const_ref<FobsValueType> const& fobs,
        af::const_ref<FcalcValueType> const& fcalc,
        bool compute_derivatives = false)
      {
        init(fobs, af::const_ref<WeightValueType>(0,0),
             fcalc, compute_derivatives);
      }

      FobsValueType
      correlation() const { return correlation_; }

      FobsValueType
      target() const { return target_; }

      af::shared<FcalcValueType>
      derivatives() const { return derivatives_; }

    protected:
      FobsValueType correlation_;
      FobsValueType target_;
      af::shared<FcalcValueType> derivatives_;

      void init(
        af::const_ref<FobsValueType> const& fobs,
        af::const_ref<WeightValueType> const& weights,
        af::const_ref<FcalcValueType> const& fcalc,
        bool compute_derivatives);
  };

  template <typename FobsValueType,
            typename WeightValueType,
            typename FcalcValueType,
            typename SumWeightsType>
  void
  intensity_correlation<FobsValueType,
                        WeightValueType,
                        FcalcValueType,
                        SumWeightsType>
  ::init(
    af::const_ref<FobsValueType> const& fobs,
    af::const_ref<WeightValueType> const& weights,
    af::const_ref<FcalcValueType> const& fcalc,
    bool compute_derivatives)
  {
    CCTBX_ASSERT(fobs.size() == weights.size() || weights.size() == 0);
    CCTBX_ASSERT(fobs.size() == fcalc.size());
    SumWeightsType sum_weights(0);
    FobsValueType sum_x(0);
    FobsValueType sum_x2(0);
    FobsValueType sum_y(0);
    FobsValueType sum_y2(0);
    FobsValueType sum_xy(0);
    FobsValueType w(1);
    for(std::size_t i=0;i<fobs.size();i++) {
      if (weights.size()) {
        w = weights[i];
        sum_weights += weights[i];
      }
      FobsValueType x = fobs[i] * fobs[i];
      FobsValueType y = std::norm(fcalc[i]);
      sum_x += w * x;
      sum_x2 += w * x * x;
      sum_y += w * y;
      sum_y2 += w * y * y;
      sum_xy += w * x * y;
    }
    if (!weights.size()) sum_weights = fobs.size();
    FobsValueType sum_w(sum_weights);
    FobsValueType x2xx = sum_x2 - sum_x * sum_x / sum_w;
    FobsValueType y2yy = sum_y2 - sum_y * sum_y / sum_w;
    FobsValueType xyxy = sum_xy - sum_x * sum_y / sum_w;
    FobsValueType correlation_denom2 = x2xx * y2yy;
    correlation_ = 1;
    if (compute_derivatives) {
      derivatives_ = af::shared<FcalcValueType>(fobs.size());
    }
    if (correlation_denom2 > 0) {
      FobsValueType correlation_denom = std::sqrt(correlation_denom2);
      correlation_ = xyxy / correlation_denom;
      if (compute_derivatives) {
        FobsValueType two_w(2);
        for(std::size_t i=0;i<fobs.size();i++) {
          if (weights.size()) {
            two_w = 2 * weights[i];
          }
          FobsValueType x = fobs[i] * fobs[i];
          FobsValueType y = std::norm(fcalc[i]);
          FobsValueType factor_deriv =
              (y - sum_y / sum_w) * correlation_ / y2yy
            - (x - sum_x / sum_w) / correlation_denom;
          derivatives_[i] = std::conj(fcalc[i]) * two_w * factor_deriv;
        }
      }
    }
    target_ = 1 - correlation_;
  }

}}} // namespace cctbx::xray::targets

#endif // CCTBX_XRAY_TARGETS_H
