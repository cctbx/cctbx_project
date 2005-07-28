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

// ///////////////////////////////////////////////
inline
   double SIMILAR(double Y)
   {
      double EPSILON = 1.0e-15;
      double LOWERLIM = 20.0;
      double ZERO = 0.0;
      double ONE = 1.0;
      double TWO = 2.0;
      double FOUR = 4.0;
      int MAXTERMS = 150;
      double X,DPN, TOT0, TOT1, SUBTOT0, SUBTOT1;
      int N;

      X = std::abs(Y);
      TOT0 = ONE;
      SUBTOT0 = ONE;
      TOT1 = ONE;
      SUBTOT1 = ONE;
      if (X < LOWERLIM) {
          N=1;
          while ((N <= MAXTERMS) && (SUBTOT0 >= EPSILON)) {
              DPN = float(N);
              SUBTOT0 = X*X*SUBTOT0/(FOUR*DPN*DPN);
              SUBTOT1 = X*X*SUBTOT1/(FOUR*DPN*(DPN+ONE));
              TOT0 = TOT0 + SUBTOT0;
              TOT1 = TOT1 + SUBTOT1;
              N=N+1;
          }
          TOT0 = TOT1*X/(TWO*TOT0);
      }
      else {
          N=1;
          while ((N <= MAXTERMS) && (std::abs(SUBTOT0) >= EPSILON)) {
              DPN = TWO*float(N);
              SUBTOT0 = (DPN - ONE)*(DPN - ONE) / (FOUR*X*DPN)*SUBTOT0;
              TOT0 = TOT0 + SUBTOT0;
              TOT1 = (TWO/(ONE - DPN) - ONE) * SUBTOT0 + TOT1;
              N=N+1;
          }
          TOT0 = TOT1/TOT0;
      }
      if (Y < ZERO) TOT0 = -TOT0;
      return TOT0;
     }
//////////////////////////////////////////////////

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
inline
double maximum_likelihood_target_one_h(double fo,
                                       double fc,
                                       double a,
                                       double b,
                                       double k,
                                       int e,
                                       int c)
{
  CCTBX_ASSERT( (c == 1 || c == 0) && (e > 0) );
  //CCTBX_ASSERT( fo >= 0 && fc >= 0 );
  //CCTBX_ASSERT( std::abs(k) > 1.e-9 );
  if(k <= 0.0) k = 1.0;
  double target = 0.0;
  if(a <= 0.0 || b <= 1.e-3 || fo <= 0.0 || fc <= 0.0) {
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
inline
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
  CCTBX_ASSERT( fo >= 0 && fc > 0);
  //CCTBX_ASSERT( std::abs(k) > 1.e-9 );
  if(k <= 0.0) k = 1.0;
  std::complex<double> d_target_over_fc = std::complex<double> (0.0,0.0);
  if(a <= 0.0 || b <= 1.e-3) {
     return std::complex<double> (0.0,0.0);
  }
  a *= k;
  b *= k*k;
  double eb = e * b;
  if(c == 0) {
    double d1 = 2. * a * a * fc / eb;
    double d2 = -2. * a * fo / eb * scitbx::math::bessel::i1_over_i0(2.*a*fo*fc/eb);
    //double d2 = -2. * a * fo / eb * SIMILAR(2.*a*fo*fc/eb);
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
inline
double d_maximum_likelihood_target_one_h_over_k(double fo,
                                                double fc,
                                                double a,
                                                double b,
                                                double k,
                                                int e,
                                                int c)
{
  CCTBX_ASSERT( (c == 1 || c == 0) && (e > 0) );
  //CCTBX_ASSERT( fo >= 0 );
  //CCTBX_ASSERT( fc >  0 );
  //CCTBX_ASSERT( std::abs(k) > 1.e-9 );
  double d_target_over_k = 0.0;
  if(a <= 0.0 || b <= 1.e-10 || fo <= 0.0 || fc <= 0.0 || std::abs(k) < 1.e-10) {
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
            typename CentricValueType = bool>
  class maximum_likelihood_criterion {

  public:
    maximum_likelihood_criterion(af::const_ref<FobsValueType> const& fobs,
                                 af::const_ref<FcalcValueType> const& fcalc,
                                 af::const_ref<FobsValueType> const& alpha,
                                 af::const_ref<FobsValueType> const& beta,
                                 FobsValueType const& k,
                                 af::const_ref<EpsilonValueType> const& eps,
                                 af::const_ref<CentricValueType> const& cs,
                                 bool compute_derivatives)
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
            af::const_ref<CentricValueType> const& cs,
            bool compute_derivatives)
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
      int c = static_cast<int>(cs[i]);
      target_ += maximum_likelihood_target_one_h(fo,fc,a,b,k,e,c);
      if(compute_derivatives) {
        derivatives_[i] = d_maximum_likelihood_target_one_h_over_fc(fo,fcalc[i],a,b,k,e,c) * (1./ fobs.size());
      }
    }
    target_ /= fobs.size();
  }
  };

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
  const af::tiny<double, 4>* cos_sin_table,
  int n_steps,
  double step_for_integration,
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
      target *= step_for_integration;
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
mlhl_d_target_dfcalc_one_h(
  double fo,
  double fc,
  double pc,
  double ac,
  double bc,
  double alpha,
  double beta,
  double k,
  int epsilon,
  bool cf,
  cctbx::hendrickson_lattman<double> const& abcd,
  const af::tiny<double, 4>* cos_sin_table,
  int n_steps,
  double step_for_integration,
  std::complex<double> const& fc_complex,
  double* workspace)
{
  const double small = 1.e-9;
  if (fc < small || alpha <= 0 || beta <= 0) {
    return std::complex<double>(0,0);
  }
  alpha *= k;
  beta *= k*k;
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
        double sim = SIMILAR(val);
        derfc = sim*arg*(arg*fc + hl_a*cos_pc + hl_b*sin_pc)/val;
        derpc = sim*arg*fc*(hl_a*sin_pc - hl_b*cos_pc)/val;
      }
    }
    // calculate derivative numerically
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
      target *= step_for_integration;
      target = -std::log(target) - maxv;
      double deranot = 0;
      double derbnot = 0;
      for(int i=0;i<n_steps;i++) {
        double exp_t_w = std::exp(target+workspace[i]);
        const double* tab = cos_sin_table[i].begin();
        deranot += tab[0] * exp_t_w;
        derbnot += tab[1] * exp_t_w;
      }
      deranot *= step_for_integration;
      derbnot *= step_for_integration;
      derfc = arg*(deranot*cos_pc + derbnot*sin_pc);
      derpc = arg*(deranot*sin_pc - derbnot*cos_pc)*fc;
    }
    derfc = 2*alpha*alpha*fc/(beta*epsilon) - derfc;
  }
  // centric reflection
  else {
    double var = beta*epsilon;
    double arg = hl_a*cos_pc + hl_b*sin_pc + fo*alpha*fc/var;
    double exp_2_arg = std::exp(-2*arg);
    double tmp_tanh = (1-exp_2_arg) / (1+exp_2_arg);
    derfc = alpha*alpha*fc/var - tmp_tanh*fo*alpha/var;
    derpc = 2*tmp_tanh*(hl_a*sin_pc - hl_b*cos_pc);
  }
  return std::complex<double>(
     (derfc*ac - derpc*bc/fc)/fc,
    -(derfc*bc + derpc*ac/fc)/fc);
}

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
            typename CentricValueType = bool,
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

    af::shared<FobsValueType>
      targets() const { return targets_; }

  protected:
    FobsValueType target_;
    af::shared<FcalcValueType> derivatives_;
    af::shared<FobsValueType> targets_;
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
    targets_ = af::shared<FobsValueType>(fobs.size());
    if (compute_derivatives) {
      derivatives_ = af::shared<FcalcValueType>(fobs.size());
    }
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
    std::vector<double> workspace(n_steps);
    for(std::size_t i_h=0;i_h<fobs.size();i_h++) {
      double fo = fobs[i_h];
      double fc = std::abs(fcalc[i_h]);
      double pc = std::arg(fcalc[i_h]);
      double ac = std::real(fcalc[i_h]);
      double bc = std::imag(fcalc[i_h]);
      double tmp1 = mlhl_target_one_h(
        fo,
        fc,
        pc,
        alpha[i_h],
        beta[i_h],
        1.0,
        eps[i_h],
        cs[i_h],
        abcd[i_h],
        &*cos_sin_table.begin(),
        n_steps,
        step_for_integration,
        &*workspace.begin());
      target_ += tmp1;
      targets_[i_h] = tmp1;
      if(compute_derivatives) {
        derivatives_[i_h] = mlhl_d_target_dfcalc_one_h(
          fo,
          fc,
          pc,
          ac,
          bc,
          alpha[i_h],
          beta[i_h],
          1.0,
          eps[i_h],
          cs[i_h],
          abcd[i_h],
          &*cos_sin_table.begin(),
          n_steps,
          step_for_integration,
          fcalc[i_h],
          &*workspace.begin());
      }
    } // end loop over reflections
  }
  };
// max-like_hl end

//! Least-squares target and its derivatives w.r.t. fc
//! ls = sum( w*(fo - k*fc)**2 ) / sum( w*fo**2 )
//! scale k can be fixed or recalculated
class ls_target_with_scale_k1 {
public:
   ls_target_with_scale_k1(af::const_ref<double> const& fo,
                           af::const_ref<double> const& w,
                           af::const_ref< std::complex<double> > const& fc,
                           bool const& compute_derivatives,
                           bool const& fix_scale,
                           double const& scale=0.0)

   {
       CCTBX_ASSERT(fo.size() == fc.size() && fo.size() == w.size());
       double num         = 0.0;
       double denum       = 0.0;
       double sum_w_fo_sq = 0.0;
       for(std::size_t i=0; i < fo.size(); i++) {
           double fc_abs = std::abs(fc[i]);
           num += fo[i] * fc_abs * w[i];
           denum += fc_abs * fc_abs * w[i];
           sum_w_fo_sq += w[i]*fo[i]*fo[i];
       }
       if(fix_scale) {
          CCTBX_ASSERT(scale > 0.0);
          scale_ = scale;
       }
       else {
          CCTBX_ASSERT(scale == 0.0);
          CCTBX_ASSERT(denum > 0.0);
          scale_ = num/denum;
       }
       CCTBX_ASSERT(sum_w_fo_sq > 0.0);
       if(compute_derivatives) {
          derivatives_ = af::shared<std::complex<double> >(fo.size());
       }
       target_ = 0.0;
       for(std::size_t i=0; i < fo.size(); i++) {
          double fc_abs = std::abs(fc[i]);
          double delta = fo[i] - scale_ * fc_abs;
          target_ += w[i] * delta * delta;
          if(compute_derivatives && fc_abs != 0) {
             derivatives_[i] = -2. * scale_ * w[i] * delta
                        / (sum_w_fo_sq * fc_abs) * std::conj(fc[i]);
          }
       }
       target_ /= sum_w_fo_sq;
   }

   double target() const { return target_; }
   af::shared<std::complex<double> > derivatives() { return derivatives_; }
   double scale() const { return scale_; }
private:
   double target_, scale_;
   af::shared<std::complex<double> > derivatives_;
};

//! Least-squares target and its derivatives w.r.t. fc
//! ls = sum( w*(k*fo - fc)**2 ) / sum( w*fo**2 )
//! scale k can be fixed or recalculated
class ls_target_with_scale_k2 {
public:
   ls_target_with_scale_k2(af::const_ref<double> const& fo,
                           af::const_ref<double> const& w,
                           af::const_ref< std::complex<double> > const& fc,
                           bool const& compute_derivatives,
                           bool const& fix_scale,
                           double const& scale=0.0)

   {
       CCTBX_ASSERT(fo.size() == fc.size() && fo.size() == w.size());
       double num         = 0.0;
       double denum       = 0.0;
       double sum_w_fo_sq = 0.0;
       for(std::size_t i=0; i < fo.size(); i++) {
           double fc_abs = std::abs(fc[i]);
           num += fo[i] * fc_abs * w[i];
           denum += fo[i] * fo[i] * w[i];
           sum_w_fo_sq += w[i] * fo[i] * fo[i];
       }
       if(fix_scale) {
          CCTBX_ASSERT(scale > 0.0);
          scale_ = scale;
       }
       else {
          CCTBX_ASSERT(scale == 0.0);
          CCTBX_ASSERT(denum > 0.0);
          scale_ = num/denum;
       }
       CCTBX_ASSERT(sum_w_fo_sq > 0.0);
       if(compute_derivatives) {
          derivatives_ = af::shared<std::complex<double> >(fo.size());
       }
       target_ = 0.0;
       for(std::size_t i=0; i < fo.size(); i++) {
          double fc_abs = std::abs(fc[i]);
          double delta = scale_ * fo[i] - fc_abs;
          target_ += w[i] * delta * delta;
          if(compute_derivatives && fc_abs != 0) {
             derivatives_[i] = -2. * w[i] * delta
                        / (sum_w_fo_sq * fc_abs) * std::conj(fc[i]);
          }
       }
       target_ /= sum_w_fo_sq;
   }

   double target() const { return target_; }
   af::shared<std::complex<double> > derivatives() { return derivatives_; }
   double scale() const { return scale_; }
private:
   double target_, scale_;
   af::shared<std::complex<double> > derivatives_;
};

}}} // namespace cctbx::xray::targets

#endif // CCTBX_XRAY_TARGETS_H
