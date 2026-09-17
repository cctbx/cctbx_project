#pragma once
/// Weighting schemes for L.S. (least_squares.h)

#include <scitbx/array_family/shared.h>
#include <cctbx/miller.h>
#include <cctbx/uctbx.h>
#include <smtbx/error.h>
#include <smtbx/import_scitbx_af.h>

#include <algorithm>
#include <cmath>

namespace smtbx { namespace refinement { namespace least_squares {
  using namespace cctbx;
  using namespace cctbx::uctbx;

  template<typename T>
  class IWeightingScheme {
  public:
    virtual ~IWeightingScheme() {}
    virtual T compute(T fo_sq, T sigma, T fc_sq,
      miller::index<> const &h, boost::optional<T> scale_factor) const = 0;

    T operator()(T fo_sq, T sigma, T fc_sq,
      miller::index<> const& h, boost::optional<T> scale_factor) const
    {
      return compute(fo_sq, sigma, fc_sq, h, scale_factor);
    }
  };
  /// Given the list of computed and measured structure factors, with
  /// their e.s.d., and the scale factor K (fo_sq ~ K fc_sq),
  /// return the list of weights.
  template<typename T>
  af::shared<T> weights(IWeightingScheme<T> const &weighting_scheme,
                        af::const_ref<T> const &fo_sq,
                        af::const_ref<T> const &sigmas,
                        af::const_ref<T> const &fc_sq,
                        af::const_ref<miller::index<> > const& indices,
                        T scale_factor)
  {
    af::shared<T> result(fo_sq.size());
    for (int i=0; i<fo_sq.size(); i++) {
      result[i] = weighting_scheme(fo_sq[i], sigmas[i], fc_sq[i],
        indices[i], scale_factor);
    }
    return result;
  }

  /// Pureley statistical weighting scheme
  /** i.e. The weight for a datum is \f$1/\sigma^2\f$ */
  // Similar to W5 in https://doi.org/10.1107/S0021889899007177 with a = 1.0
  template <typename T>
  struct sigma_weighting : public IWeightingScheme<T> {
    T compute(T fo_sq, T sigma, T fc_sq,
      miller::index<> const& h, boost::optional<T> scale_factor) const
    {
      SMTBX_ASSERT(sigma > 0);
      return std::pow(sigma, -2);
    }
  };

  /// Weights uniformly equal to 1
  // Equivalent to W6 in https://doi.org/10.1107/S0021889899007177 with a = 1.0
  template <typename T>
  struct unit_weighting : public IWeightingScheme<T> {
    T compute(T fo_sq, T sigma, T fc_sq,
      miller::index<> const& h, boost::optional<T> scale_factor) const
    {
      return 1.;
    }
  };

  /// Weights based on the resolution in sin(theta/lambda)
  // Nr W3 from https://doi.org/10.1107/S0021889899007177
  template <typename T>
  struct stl_weighting : public IWeightingScheme<T> {
    T a;
    stl_weighting(unit_cell const& uc, T a = 1.0)
      : uc(uc), a(a)
    {}

    T compute(T fo_sq, T sigma, T fc_sq,
      miller::index<> const& h,
      boost::optional<T> scale_factor) const
    {
      T val = std::sqrt(uc.d_star_sq(h) / 4);
      return std::pow(val, a);
    }
    unit_cell const& uc;
  };

  /// Weighting scheme used in George Sheldrick's SHELXL
  /** The weights are function of \f$F_o^2, F_c^2, \sigma(F_o^2)\f$
      and the scale factor K:

      \f[ w = \frac{ 1 }{ \sigma(F_o^2) + (a P)^2 + K b P) } \f]

      where

      \f[ P = 1/3 \max(F_o^2, 0) + 2/3 K F_c^2 \f]

    The formula differs from that in SHELXL documentation because SHELXL
    uses \f$ F_o^2/K, \sigma(F_o^2)/K and F_c^2 \f$ in its full computation
    of the weighted square of the residual:

      \f[ w_{SHELXL}(F_o^2/K, \sigma/K, F_c^2) (F_o^2/K - F_c^2)^2
         = w (F_o^2 - K F_c^2)^2 \f]
   */
  template <typename T>
  struct mainstream_shelx_weighting : public IWeightingScheme<T> {
    T a, b;

    mainstream_shelx_weighting(T a=0.1, T b=0)
      : a(a), b(b)
    {}

    T compute(T fo_sq, T sigma, T fc_sq,
      miller::index<> const& h, boost::optional<T> scale_factor) const
    {
      SMTBX_ASSERT(scale_factor);
      T k = *scale_factor;
      T p = (std::max(fo_sq, 0.) + 2 * k * fc_sq) / 3.;
      return 1. / (sigma * sigma + std::pow(a * p, 2) + b * k * p);
    }
  };

  /// Weighting scheme used in George Sheldrick's SHELXL
  /** The weights are function of \f$F_o^2, F_c^2, \sigma(F_o^2)\f$
      and the scale factor K:

      \f[ w = \frac{ 1 }{ \sigma(F_o^2) + (a P)^2 + K b P) } \f]

      where

      \f[ P = 1/3 \max(F_o^2, 0) + 2/3 K F_c^2 \f]

    The formula differs from that in SHELXL documentation because SHELXL
    uses \f$ F_o^2/K, \sigma(F_o^2)/K and F_c^2 \f$ in its full computation
    of the weighted square of the residual:

      \f[ w_{SHELXL}(F_o^2/K, \sigma/K, F_c^2) (F_o^2/K - F_c^2)^2
         = w (F_o^2 - K F_c^2)^2 \f]
   */
  template <typename T>
  struct new_shelx_weighting : public IWeightingScheme<T> {
    T a, b;

    new_shelx_weighting(T a=0.1, T b=0)
      : a(a), b(b)
    {}

    T compute(T fo_sq, T sigma, T fc_sq,
      miller::index<> const& h, boost::optional<T> scale_factor) const
    {
      SMTBX_ASSERT(scale_factor);
      T k = *scale_factor;
      T p = (std::max(fo_sq, 0.) + 2 * k * fc_sq) / 3.;
      return 1. / (sigma * sigma + std::pow(a * p, 2) + b * k * p);
    }
  };


}}}
