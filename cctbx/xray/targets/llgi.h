#ifndef CCTBX_XRAY_TARGETS_LLGI_H
#define CCTBX_XRAY_TARGETS_LLGI_H

#include <scitbx/constants.h>
#include <scitbx/math/utils.h>

namespace cctbx { namespace xray { namespace targets { namespace llgi {

  //! Placeholder for Iobs-based maximum-likelihood target function and
  //! gradients, using experimental sigmas.
  //! Code below is copy-paste from mlf.h

  inline
  double
  target_one_h(
    double Fe,     //"F_effective"
    double Do,     //"D_Obs"
    double ESN,    //"epsilon * Sigma_N"
    double Fc,     //"F_Calc"
    double siga,   //"sigma_A"
    double RSNoSP, //"root Sigma_N over Sigma_P"
    bool centric)
  {
    double target = 0.0;

    if (!centric) {
      //acentric:
      //p_a(Fe|Fc) =
      //            2*Fe / V
      //          * exp( - Fe^2 / V )
      //          * exp( - (Do*siga*RSNoSP*Fc)^2 / V)
      //          * I_0( 2*(Fe)*(Do*siga*RSNoSP*Fc) / V )
      //
      //           = 2*Fe/V * exp(-Fe^2 /V) * exp(-A^2 / V) * I_0(2*A / V)
      // where
      // V = ESN * ( 1 - (Do * siga)^2 ) "variance"
      // A = (Do*siga*RSNoSP*Fc)^2
      //
      // -log(p(Fe|Fc)) = -log(2*Fe/V) + Fe^2/V + A^2/V - log(I_0(2*A / V))

      double A = Do * siga * RSNoSP * Fc;
      double V = ESN * ( 1 - (Do*Do * siga*siga) );
      CCTBX_ASSERT(V > 0);

      target = -std::log(2*Fe/V) + (Fe*Fe/V) + (A*A/V) + -scitbx::math::bessel::ln_of_i0(2.0*A / V);
      //target = -std::log(A/V) + (B/V) + fast_alogchI0.get(C/V) //phasertng fast method
    }
    else {
      //centric
      //p_c(Fe|Fc) =
      //            2*Fe / V
      //          * exp( - Fe^2 / V )
      //          * exp( - (Do*siga*RSNoSP*Fc)^2 / V)
      //          * I_0( 2*(Fe)*(Do*siga*RSNoSP*Fc) / V )
      //
      //           = 2*Fe/V * exp(-Fe^2 /V) * exp(-A^2 / V) * I_0(2*A / V)
      // where
      // V = ESN * ( 1 - (Do * siga)^2 ) "variance"
      // A = (Do*siga*RSNoSP*Fc)^2
      //
      // -log(p(Fe|Fc)) = -log(2*Fe/V) + Fe^2/V + A^2/V - log(I_0(2*A / V))
      double A = Do * siga * RSNoSP * Fc;
      double V = ESN * ( 1 - (Do*Do * siga*siga) );
      CCTBX_ASSERT(V > 0);

      target = -std::log(2*Fe/V) + (Fe*Fe/V) + (A*A/V) + -scitbx::math::bessel::ln_of_i0(2.0*A / V);

    }

    return target;
  }

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

  class target_and_gradients : public common_results
  {
    public:
      target_and_gradients(
        af::const_ref<double> const& f_obs,
        af::const_ref<bool> const& r_free_flags,
        af::const_ref<std::complex<double> > const& f_calc,
        af::const_ref<double> const& llgi_sigmaa,
        double scale_factor,
        af::const_ref<double> const& epsilons,
        af::const_ref<bool> const& centric_flags,
        bool compute_gradients)
      :
        common_results(f_obs.size())
      {
        std::vector<double> alpha(f_obs.size());
        std::vector<double> beta(f_obs.size());
        throw std::runtime_error("llgi target is not implemented (yet)!");
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

}}}} // namespace cctbx::xray::targets::llgi

#endif // GUARD
