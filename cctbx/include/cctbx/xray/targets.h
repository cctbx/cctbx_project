#ifndef CCTBX_XRAY_TARGETS_H
#define CCTBX_XRAY_TARGETS_H

#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <complex>
#include <cmath>
#include <scitbx/math/bessel.h>
#include <cctbx/hendrickson_lattman.h>
#include <boost/scoped_array.hpp>
#include <boost/optional.hpp>

namespace cctbx { namespace xray {

/// X-ray target function of structure factors namespace
namespace targets {

  class ls_with_scale
  {
    public:
      ls_with_scale(
        bool apply_scale_to_f_calc,
        af::const_ref<double> const& f_obs,
        af::const_ref<double> const& weights,
        af::const_ref<bool> const& r_free_flags,
        af::const_ref<std::complex<double> > const& f_calc,
        bool compute_gradients,
        double scale_factor)
      :
        apply_scale_to_f_calc_(apply_scale_to_f_calc)
      {
        CCTBX_ASSERT(weights.size() == f_obs.size());
        CCTBX_ASSERT(r_free_flags.size() == 0
                  || r_free_flags.size() == f_obs.size());
        CCTBX_ASSERT(f_calc.size() == f_obs.size());
        CCTBX_ASSERT(scale_factor >= 0);
        const bool* rff = r_free_flags.begin();
        double num = 0;
        double denom = 0;
        std::size_t n_work = 0;
        double sum_w_fo_sq_work = 0;
        double sum_w_fo_sq_test = 0;
        // scale factor computation uses all f_obs and f_calc
        for(std::size_t i=0;i<f_obs.size();i++) {
          double fc_abs = std::abs(f_calc[i]);
          double w_fo = weights[i] * f_obs[i];
          num += w_fo * fc_abs;
          double w_fo_sq = w_fo * f_obs[i];
          if (apply_scale_to_f_calc_) {
            denom += weights[i] * fc_abs * fc_abs;
          }
          else {
            denom += w_fo_sq;
          }
          if (!rff || !rff[i]) {
            n_work++;
            sum_w_fo_sq_work += w_fo_sq;
          }
          else {
            sum_w_fo_sq_test += w_fo_sq;
          }
        }
        if (scale_factor == 0) {
          CCTBX_ASSERT(denom > 0);
          scale_factor_ = num / denom;
        }
        else {
          scale_factor_ = scale_factor;
        }
        CCTBX_ASSERT(sum_w_fo_sq_work > 0);
        if (compute_gradients) {
          gradients_work_.reserve(n_work);
        }
        target_work_ = 0;
        double target_test = 0;
        double grad_factor = -2;
        if (apply_scale_to_f_calc_) grad_factor *= scale_factor_;
        for(std::size_t i=0;i<f_obs.size();i++) {
          double fc_abs = std::abs(f_calc[i]);
          double delta;
          if (apply_scale_to_f_calc_) {
            delta = f_obs[i] - scale_factor_ * fc_abs;
          }
          else {
            delta = scale_factor_ * f_obs[i] - fc_abs;
          }
          double wd =  weights[i] * delta;
          if (!rff || !rff[i]) {
            target_work_ += wd * delta;
            if (compute_gradients) {
              if (fc_abs == 0) {
                gradients_work_.push_back(std::complex<double>(0,0));
              }
              else {
                gradients_work_.push_back(
                  grad_factor * wd / (sum_w_fo_sq_work * fc_abs) * f_calc[i]);
              }
            }
          }
          else {
            target_test += wd * delta;
          }
        }
        target_work_ /= sum_w_fo_sq_work;
        if (rff) {
          CCTBX_ASSERT(sum_w_fo_sq_test > 0);
          target_test_ = boost::optional<double>(
            target_test / sum_w_fo_sq_test);
        }
      }

      bool
      apply_scale_to_f_calc() const { return apply_scale_to_f_calc_; }

      double
      scale_factor() const { return scale_factor_; }

      double
      target_work() const { return target_work_; }

      boost::optional<double>
      target_test() const { return target_test_; }

      af::shared<std::complex<double> > const&
      gradients_work() { return gradients_work_; }

    private:
      bool apply_scale_to_f_calc_;
      double target_work_;
      boost::optional<double> target_test_;
      double scale_factor_;
      af::shared<std::complex<double> > gradients_work_;
  };

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

  /// A functor representing a least-squares residual.
  /**
  The least-square residual is defined as
  \f[
    \frac
    {\sum_i w_i \left( F_{o,i} - k \left| F_{c,i} \right| \right)^2}
    {\sum_i w_i F_{o,i}^2}
  \f]
   where \f$F_{o,i}\f$ is the i-th observed F whereas \f$F_{c,i}\f$ is the calculated F
   corresponding to \f$F_{o,i}\f$.
   It also features the weights \f$\{w_i\}\f$ and the scale factor \f$k\f$.
  */
  template <typename FobsValueType = double,
            typename WeightValueType = FobsValueType,
            typename FcalcValueType = std::complex<FobsValueType> >
  class least_squares_residual
  {
    public:
      /// Construct an uninitialised object.
      least_squares_residual() {}

      /// Construct a weighted least-squares residual.
      /**
      @param f_obs  a reference to the array containing the observed F's
      @param weights  a reference to the array containing the weights
      @param f_calc  a reference to the array containing the calculated F's
      @param compute_derivative  whether to compute the derivatives of the residual
      w.r.t. the \f$F_{c,i}\f$
      @param scale_factor  the scale factor k; if 0 then k is computed

      */
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

      /// Construct an unit weights least-square residual.
      /** Same as the other constructor but with all weights equal to unity */
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

      /// The scale factor
      FobsValueType
      scale_factor() const { return scale_factor_; }

      /// The value of the residual
      FobsValueType
      target() const { return target_; }

      /// The vector of derivatives
      /** with respect to of the residual w.r.t. \f$F_{c,i}\f$ and only if the object
      was constructed with the flag compute_derivatives==true
      */
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
                        / (sum_w_fobs2 * abs_fcalc) * fcalc[i];
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
          derivatives_[i] = fcalc[i] * two_w * factor_deriv;
        }
      }
    }
    target_ = 1 - correlation_;
  }

  namespace detail {

    inline
    double
    similar(double y)
    {
      double epsilon = 1.0e-15;
      double lowerlim = 20.0;
      int maxterms = 150;
      double x = std::abs(y);
      double tot0 = 1;
      double subtot0 = 1;
      double tot1 = 1;
      double subtot1 = 1;
      if (x < lowerlim) {
        int n = 1;
        while ((n <= maxterms) && (subtot0 >= epsilon)) {
          double dpn = static_cast<double>(n);
          subtot0 = x*x*subtot0/(4*dpn*dpn);
          subtot1 = x*x*subtot1/(4*dpn*(dpn+1));
          tot0 += subtot0;
          tot1 += subtot1;
          n++;
        }
        tot0 = tot1*x/(2*tot0);
      }
      else {
        int n = 1;
        while ((n <= maxterms) && (std::abs(subtot0) >= epsilon)) {
          double dpn = static_cast<double>(2*n);
          subtot0 = (dpn - 1)*(dpn - 1) / (4*x*dpn)*subtot0;
          tot0 += subtot0;
          tot1 += (2/(1 - dpn) - 1) * subtot0;
          n++;
        }
        tot0 = tot1/tot0;
      }
      if (y < 0) tot0 = -tot0;
      return tot0;
    }

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
          n_work = 0;
          for(std::size_t i=0;i<n_refl;i++) {
            if (flags[i]) n_work++;
          }
          if (n_work != 0) {
            n_test = n_work - n_refl;
          }
          else {
            flags = 0;
            n_work = n_refl;
            n_test = 0;
          }
        }
      }

      bool
      is_work_refl(std::size_t i) const
      {
        if (flags) return flags[i];
        return true;
      }
    };

  } // namespace detail

  //! Amplitude based Maximum-Likelihood target for one miller index.
  /*! No phase information included.
      Pavel Afonine, 28-DEC-2004.
      fo   = |Fobs|
      fc   = |Fcalc|
      a, b = distribution parameters alpha and beta
      k    = overall scale coefficient
      e    = epsilon, statistical weight for reflection
      centric = flag (false for acentric, true for centric)
  */
  inline
  double
  maximum_likelihood_target_one_h(
    double fo,
    double fc,
    double a,
    double b,
    double k,
    int e,
    bool centric)
  {
    CCTBX_ASSERT(e > 0);
    if(k <= 0.0) k = 1.0;
    double target = 0.0;
    if (a <= 0.0 || b <= 1.e-3 || fo <= 0.0 || fc <= 0.0) {
      return 0.0;
    }
    a *= k;
    b *= k*k;
    double eb = e * b;
    if(!centric) {
      double t1 = -std::log( 2. * fo / eb );
      double t2 = fo * fo / eb;
      double t3 = (a * fc) * (a * fc) / eb;
      double t4 = -scitbx::math::bessel::ln_of_i0( 2. * a * fo * fc / eb );
      target = (t1 + t2 + t3 + t4);
    }
    else {
      double Pi = scitbx::constants::pi;
      double t1 = -0.5 * std::log(2. / (Pi * eb));
      double t2 = fo * fo / (2. * eb);
      double t3 = (a * fc) * (a * fc) / (2.0 * eb);
      double t4 = -a * fo * fc / eb
                  - std::log((1. + std::exp(-2.*a*fo*fc/eb))/2.);
      target = (t1 + t2 + t3 + t4);
    }
    return target;
  }

  /* \brief Gradient of amplitude based Maximum-Likelihood target for one
     Miller index w.r.t. Fcalc
   */
  /*! No phase information included.
      Pavel Afonine, 03-JAN-2005.
      fo   = |Fobs|
      fc   = Fcalc
      a, b = distribution parameters alpha and beta
      k    = overall scale coefficient
      e    = epsilon, statistical weight for reflection
      centric = flag (false for acentric, true for centric)
  */
  inline
  std::complex<double>
  d_maximum_likelihood_target_one_h_over_fc(
    double fo,
    std::complex<double> fc_complex,
    double a,
    double b,
    double k,
    int e,
    bool centric)
  {
    CCTBX_ASSERT(e > 0);
    CCTBX_ASSERT(fo >= 0);
    double fc = std::abs(fc_complex);
    CCTBX_ASSERT(fc > 0);
    if(k <= 0.0) k = 1.0;
    std::complex<double> d_target_over_fc(0, 0);
    if(a <= 0.0 || b <= 1.e-3) {
      return d_target_over_fc;
    }
    a *= k;
    b *= k*k;
    double eb = e * b;
    if(!centric) {
      double d1 = 2. * a * a * fc / eb;
      double d2 = -2. * a * fo / eb
                  * scitbx::math::bessel::i1_over_i0(2.*a*fo*fc/eb);
      d_target_over_fc = (d1 + d2) * ( std::conj(fc_complex) / fc );
    }
    else {
      double d1 = a * a * fc / eb;
      double d2 = -a * fo / eb * std::tanh(a * fo * fc / eb);
      d_target_over_fc = (d1 + d2) * ( std::conj(fc_complex) / fc );
    }
    return d_target_over_fc;
  }

  /*! \brief Gradient of amplitude based Maximum-Likelihood target for one
      Miller index w.r.t. k
   */
  /*! No phase information included.
      Pavel Afonine, 03-JAN-2005.
      fo   = |Fobs|
      fc   = |Fcalc|
      a, b = distribution parameters alpha and beta
      k    = overall scale coefficient
      e    = epsilon, statistical weight for reflection
      centric = flag (false for acentric, true for centric)
  */
  inline
  double
  d_maximum_likelihood_target_one_h_over_k(
    double fo,
    double fc,
    double a,
    double b,
    double k,
    int e,
    bool centric)
  {
    CCTBX_ASSERT(e > 0);
    double d_target_over_k = 0.0;
    if (   a <= 0.0 || b <= 1.e-10
        || fo <= 0.0 || fc <= 0.0
        || std::abs(k) < 1.e-10) {
       return 0.0;
    }
    double eb = e * b;
    if(!centric) {
      double d1 = 2. / k;
      double d2 = -2. * fo * fo / (eb * k*k*k);
      double d3 = 2. * a * fo * fc / (eb * k*k)
        * scitbx::math::bessel::i1_over_i0(2. * a * fo * fc / (eb * k));
      d_target_over_k = d1 + d2 + d3;
    }
    else {
      double d1 = 1. / k;
      double d2 = - fo * fo / (eb * k*k*k);
      double d3 = a * fo * fc / (eb * k*k)
        * std::tanh(a * fo * fc / (eb * k));
      d_target_over_k = d1 + d2 + d3;
    }
    return d_target_over_k;
  }

  //! maximum-likelihood target function and gradients
  /*! Pavel Afonine, 26-MAY-2004
   */
  class maximum_likelihood_criterion
  {
    public:
      maximum_likelihood_criterion(
        af::const_ref<double> const& f_obs,
        af::const_ref<bool> const& r_free_flags,
        af::const_ref<std::complex<double> > const& f_calc,
        af::const_ref<double> const& alpha,
        af::const_ref<double> const& beta,
        double scale_factor,
        af::const_ref<int> const& epsilons,
        af::const_ref<bool> const& centric_flags,
        bool compute_gradients)
      :
        target_work_(0)
      {
        CCTBX_ASSERT(r_free_flags.size() == 0
                  || r_free_flags.size() == f_obs.size());
        CCTBX_ASSERT(f_calc.size() == f_obs.size());
        CCTBX_ASSERT(alpha.size() == f_obs.size());
        CCTBX_ASSERT(beta.size() == f_obs.size());
        CCTBX_ASSERT(beta.size() == epsilons.size());
        CCTBX_ASSERT(epsilons.size() == f_obs.size());
        CCTBX_ASSERT(centric_flags.size() == f_obs.size());
        if (f_obs.size() == 0) return;
        detail::r_free_flags_stats rffs(f_obs.size(), r_free_flags.begin());
        CCTBX_ASSERT(rffs.n_work != 0);
        double one_over_n_work = 1./ rffs.n_work;
        if (compute_gradients) {
          gradients_work_.reserve(rffs.n_work);
        }
        double target_test = 0;
        for(std::size_t i=0;i<f_obs.size();i++) {
          double fo = f_obs[i];
          double fc = std::abs(f_calc[i]);
          double a  = alpha[i];
          double b  = beta[i];
          int e  = epsilons[i];
          bool c = centric_flags[i];
          double t = maximum_likelihood_target_one_h(
            fo, fc, a, b, scale_factor, e, c);
          if (rffs.is_work_refl(i)) {
            target_work_ += t;
            if (compute_gradients) {
              gradients_work_.push_back(std::conj(
                d_maximum_likelihood_target_one_h_over_fc(
                  fo, f_calc[i], a, b, scale_factor, e, c)) * one_over_n_work);
            }
          }
          else {
            target_test += t;
          }
        }
        target_work_ *= one_over_n_work;
        if (rffs.n_test != 0) {
          target_test_ = boost::optional<double>(target_test / rffs.n_test);
        }
      }

      double
      target_work() const { return target_work_; }

      boost::optional<double>
      target_test() const { return target_test_; }

      af::shared<std::complex<double> >
      gradients_work() const { return gradients_work_; }

    protected:
      double target_work_;
      boost::optional<double> target_test_;
      af::shared<std::complex<double> > gradients_work_;
  };

  namespace detail {

    inline
    double
    mlhl_target_one_h(
      double fo,
      double fc,
      double pc,
      double alpha,
      double beta,
      double k,
      int epsilon,
      bool cf,
      cctbx::hendrickson_lattman<double> const& abcd,
      const af::tiny_plain<double, 4>* cos_sin_table,
      int n_steps,
      double integration_step_size,
      double* workspace)
    {
      CCTBX_ASSERT(fo >= 0);
      CCTBX_ASSERT(fc >= 0);
      const double small = 1.e-9;
      CCTBX_ASSERT(std::abs(k) > small);
      if(alpha <= 0 || beta <= 0) return 0;
      double target = 0;
      alpha *= k;
      beta *= k*k;
      double hl_a = abcd.a();
      double hl_b = abcd.b();
      double hl_c = abcd.c();
      double hl_d = abcd.d();
      // acentric reflection
      if(!cf) {
        double arg = 2*alpha*fo*fc/(beta*epsilon);
        double a_prime = arg * std::cos(pc) + hl_a;
        double b_prime = arg * std::sin(pc) + hl_b;
        // calculate target analytically
        if(std::abs(hl_c) < small && std::abs(hl_d) < small) {
          double val = std::sqrt(a_prime*a_prime + b_prime*b_prime);
          if(std::abs(hl_a) < small && std::abs(hl_b) < small) {
            val = arg;
          }
          target = scitbx::math::bessel::ln_of_i0(val);
        }
        // calculate target numerically
        else {
          double maxv = 0;
          for(int i=0;i<n_steps;i++) {
            const double* tab = cos_sin_table[i].begin();
            double term = a_prime * tab[0]
                        + b_prime * tab[1]
                        + hl_c    * tab[2]
                        + hl_d    * tab[3];
            if (maxv < term) maxv = term;
            workspace[i] = term;
          }
          for(int i=0;i<n_steps;i++) {
            target += std::exp(-maxv+workspace[i]);
          }
          target *= integration_step_size;
          target = std::log(target) + maxv;
        }
        target = (fo*fo+alpha*alpha*fc*fc)/(beta*epsilon) - target;
      }
      // centric reflection
      else {
        double var = beta*epsilon;
        double arg = fo*alpha*fc/var;
        arg += hl_a*std::cos(pc) + hl_b*std::sin(pc);
        double mabsarg = -std::abs(arg);
        target = (fo*fo + alpha*alpha*fc*fc)/(2*var) + mabsarg
               - std::log((1 + std::exp(2*mabsarg))/2);
      }
      return target;
    }

    inline
    std::complex<double>
    mlhl_d_target_d_f_calc_one_h(
      double fo,
      double fc,
      double pc,
      double ac,
      double bc,
      double alpha,
      double beta,
      int epsilon,
      bool cf,
      cctbx::hendrickson_lattman<double> const& abcd,
      const af::tiny_plain<double, 4>* cos_sin_table,
      int n_steps,
      double integration_step_size,
      double* workspace)
    {
      const double small = 1.e-9;
      if (fc < small || alpha <= 0 || beta <= 0) {
        return std::complex<double>(0,0);
      }
      double derfc = 0;
      double derpc = 0;
      double cos_pc = std::cos(pc);
      double sin_pc = std::sin(pc);
      double hl_a = abcd.a();
      double hl_b = abcd.b();
      // acentric reflection
      if (!cf) {
        double hl_c = abcd.c();
        double hl_d = abcd.d();
        double arg = 2*alpha*fo/(beta*epsilon);
        double a_prime = arg * fc * cos_pc + hl_a;
        double b_prime = arg * fc * sin_pc + hl_b;
        if (std::abs(hl_c) < small && std::abs(hl_d) < small) {
          double val = std::sqrt(a_prime*a_prime + b_prime*b_prime);
          if(val < small) {
            derfc = 0;
            derpc = 0;
          }
          else {
            double sim = similar(val);
            derfc = sim*arg*(arg*fc + hl_a*cos_pc + hl_b*sin_pc)/val;
            derpc = sim*arg*fc*(hl_a*sin_pc - hl_b*cos_pc)/val;
          }
        }
        // calculate gradients numerically
        else {
          double maxv = 0;
          for(int i=0;i<n_steps;i++) {
            const double* tab = cos_sin_table[i].begin();
            double term = a_prime * tab[0]
                        + b_prime * tab[1]
                        + hl_c    * tab[2]
                        + hl_d    * tab[3];
            if (maxv < term) maxv = term;
            workspace[i] = term;
          }
          double target = 0;
          for(int i=0;i<n_steps;i++) {
            target += std::exp(-maxv+workspace[i]);
          }
          target *= integration_step_size;
          target = -std::log(target) - maxv;
          double deranot = 0;
          double derbnot = 0;
          for(int i=0;i<n_steps;i++) {
            double exp_t_w = std::exp(target+workspace[i]);
            const double* tab = cos_sin_table[i].begin();
            deranot += tab[0] * exp_t_w;
            derbnot += tab[1] * exp_t_w;
          }
          deranot *= integration_step_size;
          derbnot *= integration_step_size;
          derfc = arg*(deranot*cos_pc + derbnot*sin_pc);
          derpc = arg*(deranot*sin_pc - derbnot*cos_pc)*fc;
        }
        derfc = 2*alpha*alpha*fc/(beta*epsilon) - derfc;
      }
      // centric reflection
      else {
        double var = beta*epsilon;
        double arg = hl_a*cos_pc + hl_b*sin_pc + fo*alpha*fc/var;
        double mtwo_arg = -2*arg;
        if(mtwo_arg > 30.) mtwo_arg = 30.0;
        double exp_2_arg = std::exp(mtwo_arg);
        double tmp_tanh = (1-exp_2_arg) / (1+exp_2_arg);
        derfc = alpha*alpha*fc/var - tmp_tanh*fo*alpha/var;
        derpc = 2*tmp_tanh*(hl_a*sin_pc - hl_b*cos_pc);
      }
      return std::complex<double>(
         (derfc*ac - derpc*bc/fc)/fc,
        -(derfc*bc + derpc*ac/fc)/fc);
    }

  } // namespace detail

  //! Maximum-likelihood target function and gradients.
  /*! Incorporates experimental phase information as HL coefficients ABCD.
      As described by Pannu et al, Acta Cryst. (1998). D54, 1285-1294.
      All the equations are reformulated in terms of alpha/beta.
      Pavel Afonine // 14-DEC-2004
   */
  class maximum_likelihood_criterion_hl
  {
    public:
      maximum_likelihood_criterion_hl(
        af::const_ref<double> const& f_obs,
        af::const_ref<bool> const& r_free_flags,
        af::const_ref<cctbx::hendrickson_lattman<double> > const&
          experimental_phases,
        af::const_ref<std::complex<double> > const& f_calc,
        af::const_ref<double> const& alpha,
        af::const_ref<double> const& beta,
        af::const_ref<int> const& epsilons,
        af::const_ref<bool> const& centric_flags,
        double integration_step_size,
        bool compute_gradients)
      :
        target_work_(0)
      {
        CCTBX_ASSERT(r_free_flags.size() == 0
                  || r_free_flags.size() == f_obs.size());
        CCTBX_ASSERT(experimental_phases.size() == f_obs.size());
        CCTBX_ASSERT(f_calc.size() == f_obs.size());
        CCTBX_ASSERT(alpha.size() == f_obs.size());
        CCTBX_ASSERT(beta.size() == f_obs.size());
        CCTBX_ASSERT(epsilons.size() == f_obs.size());
        CCTBX_ASSERT(centric_flags.size() == f_obs.size());
        CCTBX_ASSERT(integration_step_size > 0);
        if (f_obs.size() == 0) return;
        detail::r_free_flags_stats rffs(f_obs.size(), r_free_flags.begin());
        CCTBX_ASSERT(rffs.n_work != 0);
        double one_over_n_work = 1./ rffs.n_work;
        if (compute_gradients) {
          gradients_work_.reserve(rffs.n_work);
        }
        hendrickson_lattman<double>::phase_integration_cos_sin_table
          cos_sin_table(scitbx::math::iround(360 / integration_step_size));
        CCTBX_ASSERT(cos_sin_table.n_steps > 0);
        boost::scoped_array<double> workspace(
          new double[cos_sin_table.n_steps]);
        double target_test = 0;
        for(std::size_t i=0;i<f_obs.size();i++) {
          double fo = f_obs[i];
          double fc = std::abs(f_calc[i]);
          double pc = std::arg(f_calc[i]);
          double ac = std::real(f_calc[i]);
          double bc = std::imag(f_calc[i]);
          double t = detail::mlhl_target_one_h(
            fo,
            fc,
            pc,
            alpha[i],
            beta[i],
            1.0,
            epsilons[i],
            centric_flags[i],
            experimental_phases[i],
            cos_sin_table.data.get(),
            cos_sin_table.n_steps,
            integration_step_size,
            workspace.get());
          if (rffs.is_work_refl(i)) {
            target_work_ += t;
            if(compute_gradients) {
              gradients_work_.push_back(std::conj(
                detail::mlhl_d_target_d_f_calc_one_h(
                  fo,
                  fc,
                  pc,
                  ac,
                  bc,
                  alpha[i],
                  beta[i],
                  epsilons[i],
                  centric_flags[i],
                  experimental_phases[i],
                  cos_sin_table.data.get(),
                  cos_sin_table.n_steps,
                  integration_step_size,
                  workspace.get())) * one_over_n_work);
            }
          }
          else {
            target_test += t;
          }
        }
        target_work_ *= one_over_n_work;
        if (rffs.n_test != 0) {
          target_test_ = boost::optional<double>(target_test / rffs.n_test);
        }
      }

      double
      target_work() const { return target_work_; }

      boost::optional<double>
      target_test() const { return target_test_; }

      af::shared<std::complex<double> >
      gradients_work() const { return gradients_work_; }

    protected:
      double target_work_;
      boost::optional<double> target_test_;
      af::shared<std::complex<double> > gradients_work_;
  };

}}} // namespace cctbx::xray::targets

#endif // CCTBX_XRAY_TARGETS_H
