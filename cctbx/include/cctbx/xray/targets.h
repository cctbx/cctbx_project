#ifndef CCTBX_XRAY_TARGETS_H
#define CCTBX_XRAY_TARGETS_H

#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <complex>
#include <cmath>
#include <scitbx/math/bessel.h>
#include <cctbx/hendrickson_lattman.h>


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

/*  Amplitude based Maximum-Likelihood target for one miller index.
    No phase information included.
    Pavel Afonine, 28-DEC-2004.
    fo   = |Fobs|
    fc   = |Fcalc|
    a, b = distribution parameters alpha and beta
    eps  = epsilon, statistical weight for reflection
    c    = flag (1 fro acentric, 0 for centric)
    k    = overall scale coefficient
*/
double maximum_likelihood_target_one_h(double fo,
                                       double fc,
                                       double a,
                                       double b,
                                       double k,
                                       int e,
                                       int c)
{
  CCTBX_ASSERT( (c == 1 || c == 0) && (e > 0) );
  CCTBX_ASSERT( fo >= 0 && fc >= 0 );
  CCTBX_ASSERT( std::abs(k) > 1.e-9 );
  double target = 0.0;
  if(a <= 0.0 || b <= 0.0) {
     return 0.0;
  }
  a *= k;
  b *= k*k;
  double eb = e * b;
  if(c == 0) {
    double t1 = -std::log( 2. * fo / eb );
    double t2 = fo * fo / eb;
    double t3 = (a * fc) * (a * fc) / eb;
    double t4 = -scitbx::math::bessel::ln_of_i0( 2. * a * fo * fc / eb );
    target = (t1 + t2 + t3 + t4);
  }
  if(c == 1) {
    double Pi = scitbx::constants::pi;
    double t1 = -0.5 * std::log(2. / (Pi * eb));
    double t2 = fo * fo / (2. * eb);
    double t3 = (a * fc) * (a * fc) / (2.0 * eb);
    double t4 = -a * fo * fc / eb - std::log((1. + std::exp(-2.*a*fo*fc/eb))/2.);
    target = (t1 + t2 + t3 + t4);
  }
  return target;
}

/*  Derivatiove of Amplitude based Maximum-Likelihhod target for one miller
    index w.r.t. Fcalc
    No phase information included.
    Pavel Afonine, 03-JAN-2005.
    fo   = |Fobs|
    fc   = Fcalc
    a, b = distribution parameters alpha and beta
    eps  = epsilon, statistical weight for reflection
    c    = flag (0 fro acentric, 1 for centric)
    k    = overall scale coefficient
*/
std::complex<double> d_maximum_likelihood_target_one_h_over_fc(
                                               double fo,
                                               std::complex<double> fc_complex,
                                               double a,
                                               double b,
                                               double k,
                                               int e,
                                               int c)
{
  double fc = std::abs(fc_complex);
  CCTBX_ASSERT( (c == 1 || c == 0) && (e > 0) );
  CCTBX_ASSERT( fo >= 0 && fc > 0 && k > 0. );
  std::complex<double> d_target_over_fc = std::complex<double> (0.0,0.0);
  if(a <= 0.0 || b <= 0.0) {
     return std::complex<double> (0.0,0.0);
  }
  a *= k;
  b *= k*k;
  double eb = e * b;
  if(c == 0) {
    double d1 = 2. * a * a * fc / eb;
    double d2 = -2. * a * fo / eb * scitbx::math::bessel::i1_over_i0(2.*a*fo*fc/eb);
    d_target_over_fc = (d1 + d2) * ( std::conj(fc_complex) / fc );
  }
  if(c == 1) {
    double d1 = a * a * fc / eb;
    double d2 = -a * fo / eb * std::tanh(a * fo * fc / eb);
    d_target_over_fc = (d1 + d2) * ( std::conj(fc_complex) / fc );
  }
  return d_target_over_fc;
}

/*  Derivative of Amplitude based Maximum-Likelihood target for one miller
    index w.r.t. k
    No phase information included.
    Pavel Afonine, 03-JAN-2005.
    fo   = |Fobs|
    fc   = |Fcalc|
    a, b = distribution parameters alpha and beta
    eps  = epsilon, statistical weight for reflection
    c    = flag (0 fro acentric, 1 for centric)
    k    = overall scale coefficient
*/
double d_maximum_likelihood_target_one_h_over_k(double fo,
                                                double fc,
                                                double a,
                                                double b,
                                                double k,
                                                int e,
                                                int c)
{
  CCTBX_ASSERT( (c == 1 || c == 0) && (e > 0) );
  CCTBX_ASSERT( fo >= 0 && fc > 0 );
  CCTBX_ASSERT( std::abs(k) > 1.e-9 );
  double d_target_over_k = 0.0;
  if(a <= 0.0 || b <= 0.0) {
     return 0.0;
  }
  double eb = e * b;
  if(c == 0) {
    double d1 = 2. / k;
    double d2 = -2. * fo * fo / (eb * k*k*k);
    double d3 = 2. * a * fo * fc / (eb * k*k) *
                scitbx::math::bessel::i1_over_i0(2. * a * fo * fc / (eb * k));
    d_target_over_k = d1 + d2 + d3;
  }
  if(c == 1) {
    double d1 = 1. / k;
    double d2 = - fo * fo / (eb * k*k*k);
    double d3 = a * fo * fc / (eb * k*k) * std::tanh(a * fo * fc / (eb * k));
    d_target_over_k = d1 + d2 + d3;
  }
  return d_target_over_k;
}

// maximum-likelihood target function and gradients
// Pavel Afonine, 26-MAY-2004
  template <typename FobsValueType    = double,
            typename FcalcValueType   = std::complex<FobsValueType>,
            typename EpsilonValueType = int,
            typename FlagValueType    = bool>
  class maximum_likelihood_criterion {

  public:
    maximum_likelihood_criterion(af::const_ref<FobsValueType> const& fobs,
                                 af::const_ref<FcalcValueType> const& fcalc,
                                 af::const_ref<FobsValueType> const& alpha,
                                 af::const_ref<FobsValueType> const& beta,
                                 FobsValueType const& k,
                                 af::const_ref<EpsilonValueType> const& eps,
                                 af::const_ref<EpsilonValueType> const& cs,
                                 FlagValueType const& compute_derivatives)
    {
       compute(fobs,fcalc,alpha,beta,k,eps,cs,compute_derivatives);
    }

    FobsValueType
      target() const { return target_; }

    af::shared<FcalcValueType>
      derivatives() const { return derivatives_; }

  protected:
    FobsValueType target_;
    af::shared<FcalcValueType> derivatives_;
    void
    compute(af::const_ref<FobsValueType> const& fobs,
            af::const_ref<FcalcValueType> const& fcalc,
            af::const_ref<FobsValueType> const& alpha,
            af::const_ref<FobsValueType> const& beta,
            FobsValueType const& k,
            af::const_ref<EpsilonValueType> const& eps,
            af::const_ref<EpsilonValueType> const& cs,
            FlagValueType const& compute_derivatives)
  {
    CCTBX_ASSERT(fobs.size()==fcalc.size()&&alpha.size()==beta.size());
    CCTBX_ASSERT(beta.size()==eps.size()&&eps.size()==cs.size());
    CCTBX_ASSERT(fobs.size()==alpha.size());
    target_ = 0;
    if (compute_derivatives) {
      derivatives_ = af::shared<FcalcValueType>(fobs.size());
    }
    for(std::size_t i=0;i<fobs.size();i++) {
      FobsValueType    fo = fobs[i];
      FobsValueType    fc = std::abs(fcalc[i]);
      FobsValueType    a  = alpha[i];
      FobsValueType    b  = beta[i];
      EpsilonValueType e  = eps[i];
      EpsilonValueType c  = cs[i];
      target_ += maximum_likelihood_target_one_h(fo,fc,a,b,k,e,c);
      if(compute_derivatives) {
        derivatives_[i] = d_maximum_likelihood_target_one_h_over_fc(fo,fcalc[i],a,b,k,e,c);
      }
    }
  }
  };

/*
   Maximum-likelihood target function and gradients.
   Incorporates experimental phase information as HL coefficients ABCD.
   As described by Pannu et al, Acta Cryst. (1998). D54, 1285-1294.
   All the equations are reformulated in terms of alpha/beta.
   Pavel Afonine // 14-DEC-2004
*/
  template <typename FobsValueType    = double,
            typename FcalcValueType   = std::complex<FobsValueType>,
            typename AlphaValueType   = FobsValueType,
            typename BetaValueType    = FobsValueType,
            typename EpsilonValueType = int,
            typename CentricValueType = EpsilonValueType,
            typename FlagValueType    = bool,
            typename HLValueType      = cctbx::hendrickson_lattman<FobsValueType> >
  class maximum_likelihood_criterion_hl {

  public:
    maximum_likelihood_criterion_hl(
      af::const_ref<FobsValueType> const& fobs,
      af::const_ref<FcalcValueType> const& fcalc,
      af::const_ref<AlphaValueType> const& alpha,
      af::const_ref<BetaValueType> const& beta,
      af::const_ref<EpsilonValueType> const& eps,
      af::const_ref<CentricValueType> const& cs,
      FlagValueType const& compute_derivatives,
      af::const_ref<HLValueType> const& abcd,
      FobsValueType const& step_for_integration)
    {
       compute(fobs,fcalc,alpha,beta,eps,cs,compute_derivatives,abcd,
               step_for_integration);
    }

    FobsValueType
      target() const { return target_; }

    af::shared<FcalcValueType>
      derivatives() const { return derivatives_; }

  protected:
    FobsValueType target_;
    af::shared<FcalcValueType> derivatives_;
    void
    compute(af::const_ref<FobsValueType> const& fobs,
            af::const_ref<FcalcValueType> const& fcalc,
            af::const_ref<AlphaValueType> const& alpha,
            af::const_ref<BetaValueType> const& beta,
            af::const_ref<EpsilonValueType> const& eps,
            af::const_ref<CentricValueType> const& cs,
            FlagValueType const& compute_derivatives,
            af::const_ref<HLValueType> const& abcd,
            FobsValueType const& step_for_integration)
  {
    CCTBX_ASSERT(fobs.size()==fcalc.size()&&alpha.size()==beta.size());
    CCTBX_ASSERT(beta.size()==eps.size()&&eps.size()==cs.size());
    CCTBX_ASSERT(fobs.size()==alpha.size());
    CCTBX_ASSERT(step_for_integration > 0.);
    CCTBX_ASSERT(abcd.size() == fobs.size());
    target_ = 0;
    if (compute_derivatives) {
      derivatives_ = af::shared<FcalcValueType>(fobs.size());
    }
    // precomute sin cos table   NINT(x) = int(ceil(x+0.5)-(fmod(x*0.5+0.25,1.0)!=0))
    int n_steps = 360./step_for_integration;
    double angular_step = scitbx::constants::two_pi / n_steps;
    std::vector<af::tiny<double, 4> > cos_sin_table;
    cos_sin_table.reserve(n_steps);
    for(int i_step=0;i_step<n_steps;i_step++) {
      double angle = i_step * angular_step;
      cos_sin_table.push_back(af::tiny<double, 4>(std::cos(angle),
                                                  std::sin(angle),
                                                  std::cos(angle+angle),
                                                  std::sin(angle+angle)));
    }

    for(std::size_t i_h=0;i_h<fobs.size();i_h++) {
      double fo = fobs[i_h];
      double fc = std::abs(fcalc[i_h]);
      double pc = std::arg(fcalc[i_h]);
      double ac = std::real(fcalc[i_h]);
      double bc = std::imag(fcalc[i_h]);
      CCTBX_ASSERT(std::abs(pc-std::atan2(std::imag(fcalc[i_h]), std::real(fcalc[i_h]))) < 0.001);
      double a  = alpha[i_h];
      double b  = beta[i_h];
      int    e  = eps[i_h];
      // acentric: c = 0, centric: c = 1
      int    c  = cs[i_h];
      double X  = fo * fc * a / (e * b);
      double L1 = (fo * fo + a * a * fc * fc) / (2. * e * b);
      if(c == 0) {
        X  *= 2.;
        L1 *= 2.;
      }
      double tgx = 0.0;
      double small = 1.e-6;
      double A = abcd[i_h].a();
      double B = abcd[i_h].b();
      double C = abcd[i_h].c();
      double D = abcd[i_h].d();
      double A_prime = X*std::cos(pc) + A;
      double B_prime = X*std::sin(pc) + B;
      // if C = D = 0: calculate ML target explicitly by formulas:
      // centric and acentric terms accumulated at once.
      // use ln(ch(x)) = x + ln[1 + exp(-2 * x)] - ln(2) to avoid overflow
      if((std::abs(C) < small) && (std::abs(D) < small)) {
        double arg = std::sqrt(A_prime * A_prime + B_prime * B_prime);
        tgx = (1-c)*( L1 - scitbx::math::bessel::ln_of_i0(arg) ) +
                     c  *( L1 - ( arg + std::log(1.+std::exp(-2.*arg)) - std::log(2.) ) );
        if(compute_derivatives) {
          double X_d = 2. * fo * a / (e * b);
          //if(c == 0) X_d *= 2.;
          double A_prime_d = X_d * fc * std::cos(pc) + A;
          double B_prime_d = X_d * fc * std::sin(pc) + B;
          double arg_d = std::sqrt(A_prime_d * A_prime_d + B_prime_d * B_prime_d);
          if(arg_d < small) {
            double derfc = 0.;
            double derpc = 0.;
            double d1 = derfc*ac - derpc*bc/fc;
            double d2 = derfc*bc + derpc*ac/fc;
            //derivatives_[i_h] = std::complex<double> (d1,d2)/fc;
          }
          else {
            double sim = scitbx::math::bessel::i1_over_i0(arg_d);
            double derfc = sim * X_d * (X_d*fc + A*std::cos(pc) + B*std::sin(pc))/arg_d;
            double derpc = sim * X_d * fc*(A*std::sin(pc) - B*std::cos(pc))/arg_d;
            double d1 = derfc*ac - derpc*bc/fc;
            double d2 = derfc*bc + derpc*ac/fc;
//std::cout<<" "<<std::endl;
//std::cout<<sim<<" "<<X_d<<" "<<arg_d<<" "<<derfc<<" "<<derpc<<std::endl;
//std::cout<<d1<<" "<<d2<<" "<<A_prime_d<<" "<<B_prime_d<<" "<<C<<" "<<D<<" "<<c<<std::endl;
            //derivatives_[i_h] = std::complex<double> (d1,d2)/fc;
          }
        }
      }
      // otherwise do it numerically
      else {
        // acentric reflections
        if(c == 0) {
          double maxv = 0.;
          for(int i=0;i<n_steps;i++) {
            maxv = std::max(maxv, A_prime*cos_sin_table[i][0]+
                                  B_prime*cos_sin_table[i][1]+
                                  C      *cos_sin_table[i][2]+
                                  D      *cos_sin_table[i][3]);
          }
          double integral = 0.;
          for(int i=0;i<n_steps;i++) {
            integral = integral + std::exp(-maxv+A_prime*cos_sin_table[i][0]+
                                                 B_prime*cos_sin_table[i][1]+
                                                 C      *cos_sin_table[i][2]+
                                                 D      *cos_sin_table[i][3]);
          }
          integral = integral*step_for_integration;
          integral = std::log(integral) + maxv;
          tgx = (L1 - integral);
          if(compute_derivatives) {
            double X_d = 2. * fo * a / (e * b);
            double A_prime_d = X_d * fc * std::cos(pc) + A;
            double B_prime_d = X_d * fc * std::sin(pc) + B;
            double arg_d = std::sqrt(A_prime_d * A_prime_d + B_prime_d * B_prime_d);
            double maxv = 0.;
            for(int i=0;i<n_steps;i++) {
              maxv = std::max(maxv, A_prime_d*cos_sin_table[i][0]+
                                    B_prime_d*cos_sin_table[i][1]+
                                    C        *cos_sin_table[i][2]+
                                    D        *cos_sin_table[i][3]);
            }

            double integral = 0.;
            for(int i=0;i<n_steps;i++) {
              integral = integral + std::exp(-maxv+A_prime_d*cos_sin_table[i][0]+
                                                   B_prime_d*cos_sin_table[i][1]+
                                                   C        *cos_sin_table[i][2]+
                                                   D        *cos_sin_table[i][3]);
            }
            integral = integral*step_for_integration;
            integral = -std::log(integral) - maxv;

            double deranot = 0.;
            double derbnot = 0.;
            for(int i=0;i<n_steps;i++) {
              double tmp = std::exp(-integral+A_prime_d*cos_sin_table[i][0]+
                                              B_prime_d*cos_sin_table[i][1]+
                                              C        *cos_sin_table[i][2]+
                                              D        *cos_sin_table[i][3]);
              deranot = deranot + cos_sin_table[i][0]*tmp;
              derbnot = derbnot + cos_sin_table[i][1]*tmp;
            }

            deranot = step_for_integration*deranot;
            derbnot = step_for_integration*derbnot;
            double derfc = X_d*(deranot*std::cos(pc) + derbnot*std::sin(pc));
            double derpc = X_d*(deranot*std::sin(pc) - derbnot*std::cos(pc))*fc;
            derfc = 2.*a*a*fc/(e * b) - derfc;
            double d1 = derfc*ac - derpc*bc/fc;
            double d2 = derfc*bc + derpc*ac/fc;
            //derivatives_[i_h] = std::complex<double> (d1,d2)/fc;

          }
        }
        // centric reflections
        else {
          double arg = -std::abs(A*std::cos(pc) + B*std::sin(pc) + X);
          tgx =  L1 + arg - std::log((1. + std::exp(2. * arg)) / 2.);
          if(compute_derivatives) {
            double var = e*b;
            double arg = A*std::cos(pc) + B*std::sin(pc) + fo*a*fc/var;
            double derfc = a*a*fc/var - std::tanh(arg)*fo*a/var;
            double derpc = 2.*std::tanh(arg)*(A*std::sin(pc) - B*std::cos(pc));
            double d1 = derfc*ac - derpc*bc/fc;
            double d2 = derfc*bc + derpc*ac/fc;
            //derivatives_[i_h] = std::complex<double> (d1,d2)/fc;
          }
        }
      }
    target_ += mlhl_target_one_h(fo,
                                    fc,
                                    pc,
                                    a,
                                    b,
                                    1.0,
                                    e,
                                    c,
                                    abcd[i_h],
                                    cos_sin_table,
                                    n_steps,
                                    step_for_integration);
    if (compute_derivatives) {
      derivatives_[i_h] = mlhl_d_target_dfcalc_one_h(
                                  fo,
                                  fc,
                                  pc,
                                  ac,
                                  bc,
                                  a,
                                  b,
                                  1.0,
                                  e,
                                  c,
                                  abcd[i_h],
                                  cos_sin_table,
                                  n_steps,
                                  step_for_integration);
    }
    } // end loop over reflections
  }
  };

double mlhl_target_one_h(double fo,
                         double fc,
                         double pc,
                         double alpha,
                         double beta,
                         double k,
                         int epsilon,
                         int cf,
                         cctbx::hendrickson_lattman<double> abcd,
                         std::vector<af::tiny<double, 4> > cos_sin_table,
                         int n_steps,
                         double step_for_integration)
{
  double small = 1.e-9;
  CCTBX_ASSERT( (cf == 1 || cf == 0) && (epsilon > 0) );
  CCTBX_ASSERT( fo >= 0 && fc >= 0 );
  CCTBX_ASSERT( std::abs(k) > small );
  double target = 0.0;
  if(alpha <= 0.0 || beta <= 0.0) {
     return 0.0;
  }
  alpha *= k;
  beta *= k*k;
  double A = abcd.a();
  double B = abcd.b();
  double C = abcd.c();
  double D = abcd.d();

  // acentric reflection
  if(cf == 0) {
     double arg = 2.0*alpha*fo*fc/(beta*epsilon);
     double A_prime = arg * std::cos(pc) + A;
     double B_prime = arg * std::sin(pc) + B;
     // calculate target analytically
     if((std::abs(C) < small) && (std::abs(D) < small)) {
        double val = std::sqrt(A_prime*A_prime + B_prime*B_prime);
        target = scitbx::math::bessel::ln_of_i0(val);
     }
     // calculate target numerically
     else {
       double maxv = 0.;
       for(int i=0;i<n_steps;i++) {
           maxv = std::max(maxv, A_prime*cos_sin_table[i][0]+
                                 B_prime*cos_sin_table[i][1]+
                                 C      *cos_sin_table[i][2]+
                                 D      *cos_sin_table[i][3]);
       }
       target = 0.;
       for(int i=0;i<n_steps;i++) {
           target += std::exp(-maxv+A_prime*cos_sin_table[i][0]+
                                    B_prime*cos_sin_table[i][1]+
                                    C      *cos_sin_table[i][2]+
                                    D      *cos_sin_table[i][3]);
       }
       target *= step_for_integration;
       target = std::log(target) + maxv;
     }
     target = (fo*fo+alpha*alpha*fc*fc)/(beta*epsilon) - target;
  }
  // centric reflection
  if(cf == 1) {
     double var = beta*epsilon;
     double arg = A*std::cos(pc) + B*std::sin(pc) + fo*alpha*fc/var;
     double mabsarg = -std::abs(arg);
     target = (fo*fo + alpha*alpha*fc*fc)/(2.0*var) + mabsarg -
              std::log((1.0 + std::exp(2.0*mabsarg))/2.0);
  }
  return target;
}

std::complex<double> mlhl_d_target_dfcalc_one_h(
                                  double fo,
                                  double fc,
                                  double pc,
                                  double ac,
                                  double bc,
                                  double alpha,
                                  double beta,
                                  double k,
                                  int epsilon,
                                  int cf,
                                  cctbx::hendrickson_lattman<double> abcd,
                                  std::vector<af::tiny<double, 4> > cos_sin_table,
                                  int n_steps,
                                  double step_for_integration)
{
  double small = 1.e-9;
  CCTBX_ASSERT( (cf == 1 || cf == 0) && (beta > 0.) && (epsilon > 0) );
  CCTBX_ASSERT( fo >= 0 && fc >= 0 );
  CCTBX_ASSERT( alpha >= 0. );
  CCTBX_ASSERT( std::abs(k) > small );
  std::complex<double> d_target_over_fc = std::complex<double> (0.0,0.0);
  if(alpha <= 0.0 || beta <= 0.0) {
     return std::complex<double> (0.0,0.0);
  }
  alpha *= k;
  beta *= k*k;
  double A = abcd.a();
  double B = abcd.b();
  double C = abcd.c();
  double D = abcd.d();
  double derfc = 0.0;
  double derpc = 0.0;
  // acentric reflection
  if(cf == 0) {
     double arg = 2.0*alpha*fo/(beta*epsilon);
     double A_prime = arg * std::cos(pc) + A;
     double B_prime = arg * std::sin(pc) + B;
     // calculate derivative analytically
     if((std::abs(C) < small) && (std::abs(D) < small)) {
        double val = std::sqrt(A_prime*A_prime + B_prime*B_prime);
        if(val < small) {
           derfc = 0.0;
           derpc = 0.0;
        }
        else {
           double sim = scitbx::math::bessel::i1_over_i0(val);
           derfc = sim*arg*(arg*fc + A*std::cos(pc) + B*std::sin(pc))/val;
           derpc = sim*arg*fc*(A*std::sin(pc) - B*std::cos(pc))/val;
        }
     }
     // calculate derivative numerically
     else {
       double maxv = 0.;
       for(int i=0;i<n_steps;i++) {
           maxv = std::max(maxv, A_prime*cos_sin_table[i][0]+
                                 B_prime*cos_sin_table[i][1]+
                                 C      *cos_sin_table[i][2]+
                                 D      *cos_sin_table[i][3]);
       }
       double target = 0.;
       for(int i=0;i<n_steps;i++) {
           target += std::exp(-maxv+A_prime*cos_sin_table[i][0]+
                                    B_prime*cos_sin_table[i][1]+
                                    C      *cos_sin_table[i][2]+
                                    D      *cos_sin_table[i][3]);
       }
       target *= step_for_integration;
       target = -std::log(target) - maxv;
       double deranot = 0.0;
       double derbnot = 0.0;
       for(int i=0;i<n_steps;i++) {
           double tmp = std::exp(target + A_prime*cos_sin_table[i][0]+
                                          B_prime*cos_sin_table[i][1]+
                                          C      *cos_sin_table[i][2]+
                                          D      *cos_sin_table[i][3]);
           deranot += cos_sin_table[i][0]*tmp;
           derbnot += cos_sin_table[i][1]*tmp;
       }

       deranot *= step_for_integration;
       derbnot *= step_for_integration;
       derfc = arg*(deranot*std::cos(pc) + derbnot*std::sin(pc));
       derpc = arg*(deranot*std::sin(pc) - derbnot*std::cos(pc))*fc;
     }
     derfc = 2.0*alpha*alpha*fc/(beta*epsilon) - derfc;
  }
  // centric reflection
  if(cf == 1) {
     double var = beta*epsilon;
     double arg = A*std::cos(pc) + B*std::sin(pc) + fo*alpha*fc/var;
     derfc = alpha*alpha*fc/var - std::tanh(arg)*fo*alpha/var;
     derpc = 2.0*std::tanh(arg)*(A*std::sin(pc) - B*std::cos(pc));
  }
  if(fc < small) {
     d_target_over_fc = std::complex<double> (0.0,0.0);
  }
  else {
     double d1 = derfc*ac - derpc*bc/fc;
     double d2 = derfc*bc + derpc*ac/fc;
     d_target_over_fc = std::complex<double> (d1,d2)/fc;
  }
  return std::conj(d_target_over_fc);
}

// max-like_hl end

}}} // namespace cctbx::xray::targets

#endif // CCTBX_XRAY_TARGETS_H
