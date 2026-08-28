#ifndef CCTBX_XRAY_TARGETS_LLGI_H
#define CCTBX_XRAY_TARGETS_LLGI_H

#include <scitbx/constants.h>
#include <scitbx/math/bessel.h>

namespace cctbx { namespace xray { namespace targets { namespace llgi {

  //! Log-likelihood-gain-of-intensities target for one miller index.
  /*! sigmaA-parameterised Rice-distribution likelihood-gain target, derived
      by change of variables from the E-value-scale LLGI functor in
      phasertng (phasertng::likelihood::llgi::function::target_gradient_
      hessian, phasertng/codebase/phasertng/math/likelihood/llgi/
      function.cc) to operate directly on unnormalized Feff/Fcalc, so that
      overall/anisotropic B-factor refinement (which acts on Fcalc's
      resolution/direction-dependent falloff) sees that falloff directly,
      rather than through an E-value normalization that divides it out.
      Hessian and d/dV terms present in the phasertng functor are omitted
      (see class target_and_gradients docstring below).

      Change of variables (E-scale quantity -> F-scale quantity used here):
        EOBS = FEFF/(sqrt(TEPS)*RESN)                 -> feff (unnormalized)
        Ecalc = Fcalc/sqrt(TEPS*RESN^2*ScatFrac)       -> fc   (unnormalized,
                                                                 |Fcalc|)
        alpha-equivalent d_E = Dobs*sigmaA             -> d_F = Dobs*sigmaA
                                                                 /sqrt(ScatFrac)
        V_E = TEPS - (Dobs*sigmaA)^2                   -> V_F = TEPS*RESN^2
                                                                 *V_E
      ScatFrac is the fraction of the total scattering accounted for by the
      model as a function of resolution (equivalently, the ratio of mean
      calculated to mean expected intensity in a resolution shell); it is
      passed as an explicit argument rather than folded into sigmaa, so
      that sigmaA (the model-quality parameter refined against the LLGI
      target on the test set) and ScatFrac (the resolution-dependent
      scale/anisotropy corrector, whose fitting is meant to be driven by
      this same target during overall-B refinement) stay independently
      visible and debuggable.

      feff     = effective amplitude (Feff), NOT normalized to the E-value
                 scale. Supplied by phasertng.nacelle as the FEFF column.
      fc       = |Fcalc|, NOT normalized (ordinary calculated structure
                 factor amplitude on the same scale as feff).
      dobs     = per-reflection reliability weight (Dobs), varies from
                 reflection to reflection. Supplied by phasertng.nacelle as
                 the DOBS column.
      sigmaa   = sigmaA(resolution), a smooth function of resolution only,
                 estimated within phenix.refine (see llgi_sigmaa
                 estimator), optimized against this same target on the
                 R-free set.
      scatfrac = ScatFrac(resolution), a smooth function of resolution
                 only. Passed explicitly (see above) rather than combined
                 into sigmaa.
      k        = overall scale coefficient.
      teps     = per-reflection tNCS-adjusted normalization factor (nacelle
                 TEPS column); equals 1 in the ordinary case (no tNCS).
      resn     = per-reflection Root-EpsilonSigmaN normalization (nacelle
                 RESN column), i.e. sqrt(epsilon * mean intensity at this
                 resolution).
      centric  = flag (false for acentric, true for centric).

      This target returns a log-likelihood-*gain* relative to the
      Wilson/null-hypothesis prior (function.cc's LL, which includes the
      wll = EOBS^2/teps null-hypothesis term), negated so that, like every
      other target in this file (mlf, mli, mlhl), it is a quantity for the
      refinement machinery to *minimize*.
  */
  inline
  double
  target_one_h(
    double feff,
    double fc,
    double dobs,
    double sigmaa,
    double scatfrac,
    double k,
    double teps,
    double resn,
    bool centric)
  {
    CCTBX_ASSERT(teps > 0);
    CCTBX_ASSERT(resn > 0);
    if(k <= 0.0) k = 1.0;
    if (dobs <= 0.0 || sigmaa <= 0.0 || scatfrac <= 0.0
        || feff <= 0.0 || fc <= 0.0) {
      return 0.0;
    }
    // F-scale analog of alpha (function.cc's Dobs*sigmaA, here further
    // divided by sqrt(ScatFrac) to account for the Ecalc->Fcalc rescaling).
    double d = dobs * (sigmaa / std::sqrt(scatfrac)) * k;
    double resn_sq = resn * resn;
    // v_e is the E-scale variance, function.cc's V = teps - DobsSigaSqr;
    // v is v_e rescaled to the F-scale, v = teps*resn^2*v_e.
    double v_e = teps - d * d;
    if (v_e <= 0.0) {
      // function.cc: negative_variance() guard; no contribution.
      return 0.0;
    }
    double v = teps * resn_sq * v_e;
    double teps_resn_sq = teps * resn_sq;
    double feff_sq = feff * feff;
    double ec = d * fc; // F-scale analog of EC = |DobsSigaEcalc|
    double ec_sq = ec * ec;
    double x = 2. * feff * ec / v; // function.cc: X (scale-invariant)
    // function.cc: wll = EOBS^2/teps; EOBS^2 = feff^2/(teps*resn^2), so
    // wll = feff^2/(teps^2*resn^2).
    double wll = feff_sq / (teps * teps_resn_sq);
    double ll;
    if(!centric) {
      // function.cc: log(V_E/teps) = log(v/(teps*resn^2)/teps)
      //                            = log(v/(teps^2*resn^2))
      ll = -(std::log(v / (teps * teps_resn_sq)) + (feff_sq + ec_sq) / v);
      ll += scitbx::math::bessel::ln_of_i0(x); // function.cc: tbl_alogchI0
      ll += wll;
    }
    else {
      wll /= 2.0;
      double ll_core = -(std::log(v / (teps * teps_resn_sq))
                          + (feff_sq + ec_sq) / v);
      double x_half = x / 2.0;
      // function.cc: tbl_alogch(X) approximates log(cosh(X)).
      double ln_cosh = x_half + std::log((1. + std::exp(-2. * x_half)) / 2.);
      ll = ll_core / 2.0 + ln_cosh;
      ll += wll;
    }
    // ll is the log-likelihood-gain; negate to match this file's (and
    // mlf.h's) minimize-me convention.
    return -ll;
  }

  /* \brief Gradient of the LLGI target for one Miller index w.r.t. Fcalc.
     Ported from function.cc's dLL_by_dEC (with dEC_by_dReEC = ReEC/EC,
     dEC_by_dImEC = ImEC/EC specialised to a real, positive scale term
     d = dobs*sigmaa/sqrt(scatfrac), so EC = d*fc reduces the chain rule to
     dLL/dfc = dLL/dEC * d, matching mlf.h's treatment of Fcalc as a
     complex number with |Fcalc| entering the real-valued target). Hessian
     and d/dV terms present in phasertng's functor are omitted, since
     phenix.refine's LBFGS-based target functors only consume gradients
     w.r.t. Fcalc.
  */
  inline
  std::complex<double>
  d_target_one_h_over_fc(
    double feff,
    std::complex<double> fc_complex,
    double dobs,
    double sigmaa,
    double scatfrac,
    double k,
    double teps,
    double resn,
    bool centric)
  {
    CCTBX_ASSERT(teps > 0);
    CCTBX_ASSERT(resn > 0);
    CCTBX_ASSERT(feff >= 0);
    double fc = std::abs(fc_complex);
    if(fc == 0) return std::complex<double>(0, 0);
    if(k <= 0.0) k = 1.0;
    std::complex<double> d_target_over_fc(0, 0);
    if (dobs <= 0.0 || sigmaa <= 0.0 || scatfrac <= 0.0) {
      return d_target_over_fc;
    }
    double d = dobs * (sigmaa / std::sqrt(scatfrac)) * k;
    double v_e = teps - d * d;
    if (v_e <= 0.0) {
      return d_target_over_fc;
    }
    double v = teps * resn * resn * v_e;
    double ec = d * fc; // function.cc: EC
    double x = 2. * feff * ec / v; // function.cc: X
    double dll_by_dec; // function.cc: dLL_by_dEC
    if(!centric) {
      double bess_term = scitbx::math::bessel::i1_over_i0(x);
      dll_by_dec = (2. / v) * (feff * bess_term - ec);
    }
    else {
      double bess_term = std::tanh(x / 2.0);
      dll_by_dec = (1. / v) * (feff * bess_term - ec);
    }
    double dll_by_dfc = dll_by_dec * d; // dEC_by_dfc = d
    d_target_over_fc = dll_by_dfc * ( std::conj(fc_complex) / fc );
    // Sign-flipped to match target_one_h's minimization convention.
    return -d_target_over_fc;
  }

  //! LLGI (log-likelihood-gain-of-intensities) target function and
  //! gradients.
  /*! sigmaA-parameterised alternative to the alpha/beta-based mlf target
      (cctbx/xray/targets/mlf.h), operating directly on unnormalized
      Feff/Fcalc (see target_one_h docstring for the change-of-variables
      derivation from phasertng's E-scale LLGI functor). Dobs and Feff are
      precomputed per reflection (phasertng.nacelle), along with the TEPS
      and RESN normalization factors; sigmaA and ScatFrac are smooth
      functions of resolution only, estimated/refined within phenix.refine.
  */
  class target_and_gradients : public common_results
  {
    public:
      target_and_gradients(
        af::const_ref<double> const& f_eff,
        af::const_ref<bool> const& r_free_flags,
        af::const_ref<std::complex<double> > const& f_calc,
        af::const_ref<double> const& dobs,
        af::const_ref<double> const& sigmaa,
        af::const_ref<double> const& scatfrac,
        double scale_factor,
        af::const_ref<double> const& teps,
        af::const_ref<double> const& resn,
        af::const_ref<bool> const& centric_flags,
        bool compute_gradients)
      :
        common_results(f_eff.size())
      {
        CCTBX_ASSERT(r_free_flags.size() == 0
                  || r_free_flags.size() == f_eff.size());
        CCTBX_ASSERT(f_calc.size() == f_eff.size());
        CCTBX_ASSERT(dobs.size() == f_eff.size());
        CCTBX_ASSERT(sigmaa.size() == f_eff.size());
        CCTBX_ASSERT(scatfrac.size() == f_eff.size());
        CCTBX_ASSERT(teps.size() == f_eff.size());
        CCTBX_ASSERT(resn.size() == f_eff.size());
        CCTBX_ASSERT(centric_flags.size() == f_eff.size());
        if (f_eff.size() == 0) return;
        detail::r_free_flags_stats rffs(f_eff.size(), r_free_flags.begin());
        CCTBX_ASSERT(rffs.n_work != 0);
        double one_over_n_work = 1./ rffs.n_work;
        if (compute_gradients) {
          gradients_work_.reserve(rffs.n_work);
        }
        double target_work = 0;
        double target_test = 0;
        for(std::size_t i=0;i<f_eff.size();i++) {
          double feff = f_eff[i];
          double fc = std::abs(f_calc[i]);
          double do_ = dobs[i];
          double sa = sigmaa[i];
          double sf = scatfrac[i];
          double tp = teps[i];
          double rn = resn[i];
          bool c = centric_flags[i];
          double t = target_one_h(
            feff, fc, do_, sa, sf, scale_factor, tp, rn, c);
          target_per_reflection_[i] = t;
          if (rffs.is_work_refl(i)) {
            target_work += t;
            if (compute_gradients) {
              gradients_work_.push_back(std::conj(
                d_target_one_h_over_fc(
                  feff, f_calc[i], do_, sa, sf, scale_factor, tp, rn, c))
                  * one_over_n_work);
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
