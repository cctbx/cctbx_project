#ifndef CCTBX_XRAY_TARGETS_MLF_H
#define CCTBX_XRAY_TARGETS_MLF_H

#include <scitbx/constants.h>
#include <scitbx/math/utils.h>

namespace cctbx { namespace xray { namespace targets { namespace mlf {

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
  target_one_h(
    double fo,
    double fc,
    double a,
    double b,
    double k,
    double e,
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
  d_target_one_h_over_fc(
    double fo,
    std::complex<double> fc_complex,
    double a,
    double b,
    double k,
    double e,
    bool centric)
  {
    CCTBX_ASSERT(e > 0);
    CCTBX_ASSERT(fo >= 0);
    double fc = std::abs(fc_complex);
    if(fc == 0) return std::complex<double> (0,0);
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

  //! maximum-likelihood target function and gradients
  /*! References:
      Lunin, V.Y., Afonine, P.V. & Urzhumtsev, A. (2002).
        Acta Cryst. A58, 270-282.
      Lunin, V.Y. & Skovoroda, T.P. (1995).
        Acta Cryst. A51, 880-887.
      Pavel Afonine, 26-MAY-2004
   */
  class target_and_gradients : public common_results
  {
    public:
      target_and_gradients(
        af::const_ref<double> const& f_obs,
        af::const_ref<bool> const& r_free_flags,
        af::const_ref<std::complex<double> > const& f_calc,
        af::const_ref<double> const& alpha,
        af::const_ref<double> const& beta,
        double scale_factor,
        af::const_ref<double> const& epsilons,
        af::const_ref<bool> const& centric_flags,
        bool compute_gradients)
      :
        common_results(f_obs.size())
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
        double target_work = 0;
        double target_test = 0;
        for(std::size_t i=0;i<f_obs.size();i++) {
          double fo = f_obs[i];
          double fc = std::abs(f_calc[i]);
          double a  = alpha[i];
          double b  = beta[i];
          double e  = epsilons[i];
          bool c = centric_flags[i];
          double t = target_one_h(
            fo, fc, a, b, scale_factor, e, c);
          target_per_reflection_[i] = t;
          if (rffs.is_work_refl(i)) {
            target_work += t;
            if (compute_gradients) {
              gradients_work_.push_back(std::conj(
                d_target_one_h_over_fc(
                  fo, f_calc[i], a, b, scale_factor, e, c)) * one_over_n_work);
            }
          }
          else {
            target_test += t;
          }
        }
        target_work_ = target_work * one_over_n_work;
        if (rffs.n_test != 0) {
          target_test_ = boost::optional<double>(target_test / rffs.n_test);
        }
      }
  };

}}}} // namespace cctbx::xray::targets

#endif // GUARD
