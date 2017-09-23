#ifndef CCTBX_XRAY_TARGETS_H
#define CCTBX_XRAY_TARGETS_H

#include <cctbx/xray/targets/common_results.h>
#include <cctbx/error.h>
#include <scitbx/math/bessel.h>

namespace cctbx { namespace xray { namespace targets {

  /// The modulus of complex structure factors
  /**
  This is a policy to be passed as a template argument to least_squares_residual. The type T is to be like std::complex<>
  */
  template<class T>
  struct f_calc_modulus
  {
    /// The value |f|
    static typename T::value_type value(T f) {
      return std::abs(f);
    }

    /// The derivatives of |f| wrt to respectively the real part and the imaginary part of f.
    /**
    They are packaged into the complex number of type T, respectively in its real and imaginary part.
    */
    static T derivative(T f) {
      return f / std::abs(f);
    }
  };


  /// The modulus squared of complex structure factors
  /**
    This is a policy to be passed as a template argument to least_squares_residual. The type T is to be like std::complex<>
   */
  template<class T>
  struct f_calc_modulus_square
  {
    /// The value |f|^2
    static typename T::value_type value(T f) {
      return std::norm(f);
    }

    /// The derivatives of |f|^2 wrt to respectively the real part and the imaginary part of f.
    /**
    They are packaged into the complex number of type T, respectively in its real and imaginary part.
    */
    static T derivative(T f) {
      return 2. * f;
    }
  };

  /// A functor representing a least-squares residual.
  /**
  The least-square residual is defined as
  \f[
    \frac
    {\sum_i w_i \left( y_{o,i} - k f(F_{c,i}) \right)^2}
    {\sum_i w_i y_{o,i}}
  \f]
   where \f$y_{o,i}\f$ is the i-th observed piece of data (i.e. F or F^2 in practice)
   whereas \f$F_{c,i}\f$ is the corresponding computed structure factor and \f$f\f$ is any function
   (in practice the modulus or the modulus squared), which is represented by the policy FcalcFunctor.
   It also features the weights \f$\{w_i\}\f$ and the scale factor \f$k\f$.
  */
  template <template<typename> class FcalcFunctor,
            typename YobsValueType = double,
            typename WeightValueType = YobsValueType,
            typename FcalcValueType = std::complex<YobsValueType> >
  class least_squares_residual
  {
    public:
      /// Construct an uninitialised object.
      least_squares_residual() {}

      /// Construct a weighted least-squares residual.
      /**
      @param yobs  a reference to the array containing the observed data
      @param weights  a reference to the array containing the weights
      @param fcalc  a reference to the array containing the calculated F's
      @param compute_derivative  whether to compute the derivatives of the residual
      w.r.t. the real and imaginary parts of \f$F_{c,i}\f$
      @param scale_factor  the scale factor k; if 0 then k is computed

      */
      least_squares_residual(
        af::const_ref<YobsValueType> const& yobs,
        af::const_ref<WeightValueType> const& weights,
        af::const_ref<FcalcValueType> const& fcalc,
        bool compute_derivatives=false,
        YobsValueType const& scale_factor=0)
      :
        scale_factor_(scale_factor)
      {
        init(yobs, weights, fcalc, compute_derivatives);
      }

      /// Construct an unit weights least-square residual.
      /** Same as the other constructor but with all weights equal to unity */
      least_squares_residual(
        af::const_ref<YobsValueType> const& yobs,
        af::const_ref<FcalcValueType> const& fcalc,
        bool compute_derivatives=false,
        YobsValueType const& scale_factor=0)
      :
        scale_factor_(scale_factor)
      {
        init(yobs, af::const_ref<WeightValueType>(0,0),
             fcalc, compute_derivatives);
      }

      /// The scale factor
      YobsValueType
      scale_factor() const { return scale_factor_; }

      /// The value of the residual
      YobsValueType
      target() const { return target_; }

      /// The vector of derivatives
      /**  Only if the object was constructed with the flag compute_derivatives==true. The i-th element of the returned array is a complex number whose real (resp. imaginary) part
      contains the derivative of the residual with respect to the real (resp. imaginary) part of
      of \f$F_{c,i}\f$.
      */
      af::shared<FcalcValueType>
      derivatives() const
      {
        return derivatives_;
      }

    protected:
      YobsValueType scale_factor_;
      YobsValueType target_;
      af::shared<FcalcValueType> derivatives_;

      YobsValueType
      new_scale_factor(
                       af::const_ref<YobsValueType> const& yobs,
                       af::const_ref<WeightValueType> const& weights,
                       af::const_ref<FcalcValueType> const& fcalc);

      YobsValueType
      sum_weighted_yobs_squared(
                                af::const_ref<YobsValueType> const& values,
                                af::const_ref<WeightValueType> const& weights);

      void init(
                af::const_ref<YobsValueType> const& yobs,
                af::const_ref<WeightValueType> const& weights,
                af::const_ref<FcalcValueType> const& fcalc,
                bool compute_derivatives);
  };

  template <template<typename> class FcalcFunctor,
            typename YobsValueType,
            typename WeightValueType,
            typename FcalcValueType>
  YobsValueType
  least_squares_residual<FcalcFunctor, YobsValueType, WeightValueType, FcalcValueType>
  ::new_scale_factor(
                   af::const_ref<YobsValueType> const& yobs,
                   af::const_ref<WeightValueType> const& weights,
                   af::const_ref<FcalcValueType> const& fcalc)
  {
    CCTBX_ASSERT(yobs.size() == weights.size() || weights.size() == 0);
    CCTBX_ASSERT(yobs.size() == fcalc.size());
    YobsValueType sum_w_yobs_ycalc(0);
    YobsValueType sum_w_ycalc2(0);
    WeightValueType w(1);
    for(std::size_t i=0;i<yobs.size();i++) {
      YobsValueType ycalc = FcalcFunctor<FcalcValueType>::value(fcalc[i]);
      if (weights.size()) w = weights[i];
      sum_w_yobs_ycalc += w * yobs[i] * ycalc;
      sum_w_ycalc2 += w * ycalc * ycalc;
    }
    if (sum_w_ycalc2 == 0) {
      throw cctbx::error(
                         "Cannot calculate scale factor: sum of weights * fcalc^2 == 0.");
    }
    return sum_w_yobs_ycalc / sum_w_ycalc2;
  }

  template <template<typename> class FcalcFunctor,
            typename YobsValueType,
            typename WeightValueType,
            typename FcalcValueType>
  YobsValueType
  least_squares_residual<FcalcFunctor, YobsValueType, WeightValueType, FcalcValueType>
  ::sum_weighted_yobs_squared(
                            af::const_ref<YobsValueType> const& yobs,
                            af::const_ref<WeightValueType> const& weights
                            )
  {
    CCTBX_ASSERT(yobs.size() == weights.size() || weights.size() == 0);
    YobsValueType result = 0.;
    WeightValueType w(1);
    for(std::size_t i=0;i<yobs.size();i++) {
      if (weights.size()) w = weights[i];
      result += w * yobs[i] * yobs[i];
    }
    return result;
  }


  template <template<typename> class FcalcFunctor,
            typename YobsValueType,
            typename WeightValueType,
            typename FcalcValueType>
  void
  least_squares_residual<FcalcFunctor, YobsValueType, WeightValueType, FcalcValueType>
  ::init(
    af::const_ref<YobsValueType> const& yobs,
    af::const_ref<WeightValueType> const& weights,
    af::const_ref<FcalcValueType> const& fcalc,
    bool compute_derivatives)
  {
    if (scale_factor_ == 0) {
      scale_factor_ = new_scale_factor(yobs, weights, fcalc);
    }
    YobsValueType sum_w_yobs2 = sum_weighted_yobs_squared(yobs, weights);
    if (sum_w_yobs2 == 0) {
      throw cctbx::error(
        "Cannot calculate least-squares residual:"
        " sum of weights * yobs^2 == 0.");
    }
    YobsValueType one_over_sum_w_yobs2 = 1./sum_w_yobs2;
    target_ = 0;
    if (compute_derivatives) {
      derivatives_ = af::shared<FcalcValueType>(yobs.size());
    }
    WeightValueType w(1);
    for(std::size_t i=0;i<yobs.size();i++) {
      YobsValueType ycalc = FcalcFunctor<FcalcValueType>::value(fcalc[i]);
      YobsValueType delta = yobs[i] - scale_factor_ * ycalc;
      if (weights.size()) w = weights[i];
      target_ += w * delta * delta;
      if (compute_derivatives && ycalc != 0) {
        derivatives_[i] = -2. * scale_factor_ * w * delta
                        * FcalcFunctor<FcalcValueType>::derivative(fcalc[i])
                        * one_over_sum_w_yobs2;
      }
    }
    target_ /= sum_w_yobs2;
  }

  namespace detail {

    struct r_free_flags_stats
    {
      const bool* flags;
      std::size_t n_work;
      std::size_t n_test;

      r_free_flags_stats(
        std::size_t n_refl,
        const bool* r_free_flags_begin)
      :
        flags(r_free_flags_begin)
      {
        if (!flags) {
          n_work = n_refl;
          n_test = 0;
        }
        else {
          n_test = 0;
          for(std::size_t i=0;i<n_refl;i++) {
            if (flags[i]) n_test++;
          }
          n_work = n_refl - n_test;
          if (n_test == 0) flags = 0;
        }
      }

      bool
      is_work_refl(std::size_t i) const
      {
        if (flags) return !flags[i];
        return true;
      }
    };

  } // namespace detail

  template <typename FloatType=double,
            typename ComplexType=std::complex<double> >
  class r_factor
  {
    public:
      FloatType r_;
      FloatType scale_ls_;
      FloatType scale_r_;

      r_factor() {}

      r_factor(af::const_ref<FloatType> const& fo,
               af::const_ref<ComplexType> const& fc)
      {
        CCTBX_ASSERT(fo.size()==fc.size());
        compute_scale(fo, fc);
        r_ = compute_r_factor(fo, fc, scale_r_);
      };

      FloatType compute_r_factor(af::const_ref<FloatType> const& fo,
                                 af::const_ref<ComplexType> const& fc,
                                 FloatType scale)
      {
        CCTBX_ASSERT(fo.size()==fc.size());
        FloatType num=0.0;
        FloatType denum=0.0;
        for(std::size_t i=0; i < fo.size(); i++) {
          num += std::abs(fo[i] - std::abs(fc[i]) * scale);
          denum += fo[i];
        }
        if(denum == 0) return 1.e+9;
        return num/denum;
      };

      void compute_scale(af::const_ref<FloatType> const& fo,
                         af::const_ref<ComplexType> const& fc,
                         FloatType offset_factor = 3.,
                         FloatType step_factor = 20.) {
        FloatType num=0.0;
        FloatType denum=0.0;
        for(std::size_t i=0; i < fo.size(); i++) {
          FloatType fc_abs = std::abs(fc[i]);
          num += fo[i] * fc_abs;
          denum += fc_abs * fc_abs;
        }
        scale_ls_ = (denum == 0 ? 0 : num/denum);
        FloatType scale_ = scale_ls_ - scale_ls_/offset_factor;
        FloatType r_best = compute_r_factor(fo, fc, scale_ls_);
        FloatType step = scale_ls_/step_factor;
        scale_r_ = scale_ls_;
        while(scale_ <= scale_ls_ + scale_ls_/offset_factor) {
          FloatType r_trial = compute_r_factor(fo, fc, scale_);
          if(r_trial < r_best) {
            r_best = r_trial;
            scale_r_ = scale_;
          }
          scale_ += step;
        }
      }

      FloatType value() { return r_; }
      FloatType scale_ls() { return scale_ls_; }
      FloatType scale_r() { return scale_r_; }
  };


}}} // namespace cctbx::xray::targets

#endif // CCTBX_XRAY_TARGETS_H
