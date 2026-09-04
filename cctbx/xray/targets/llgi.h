#ifndef CCTBX_XRAY_TARGETS_LLGI_H
#define CCTBX_XRAY_TARGETS_LLGI_H

#include <scitbx/constants.h>
#include <scitbx/math/bessel.h>
#include <utility>

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
    // Symmetrised form (F-scale analog of the simplification in llgi_e.h).
    // phasertng's function.cc carries the Wilson/null-hypothesis baseline
    // as a separate additive wll = EOBS^2/teps = feff^2/(teps^2*resn^2),
    // added after subtracting feff^2/V. Those two feff terms collapse:
    //
    //   -feff^2/V + feff^2/(teps^2*resn^2) = -(d^2/teps)*feff^2/V
    //
    // (using V = teps*resn^2*(teps - d^2)), so the gain follows directly
    // by scaling feff^2, mirroring how ec_sq = (d*fc)^2 scales fc^2 -- no
    // separate baseline term to add and subtract.
    //
    // IMPORTANT: on the F-scale the multiplier is d^2/TEPS, NOT the plain
    // d^2 that appears in llgi_e.h. The two coincide only at TEPS == 1
    // (which is why the E-scale form, where TEPS is identically 1, takes
    // the simpler shape). Since tNCS support (TEPS != 1) is deferred but
    // intended, using a plain d^2 here would silently bake in a wrong
    // formula in exactly the case that support is meant to enable:
    // checked numerically, at TEPS = 2 the plain-d^2 form gives -0.618
    // where the correct value is -1.395. Verified algebraically (sympy)
    // and numerically against the original three-term form over ~39k
    // randomised cases spanning TEPS in [1,3], centric and acentric,
    // agreeing to ~1e-9 relative (that residual is ln_of_i0's own
    // tabulated-approximation noise, not a difference between the forms).
    double ee_sq = (d * d / teps) * feff_sq; // symmetric partner of ec_sq
    double ll;
    if(!centric) {
      // function.cc: log(V_E/teps) = log(v/(teps*resn^2)/teps)
      //                            = log(v/(teps^2*resn^2))
      ll = -(std::log(v / (teps * teps_resn_sq)) + (ee_sq + ec_sq) / v);
      ll += scitbx::math::bessel::ln_of_i0(x); // function.cc: tbl_alogchI0
    }
    else {
      double ll_core = -(std::log(v / (teps * teps_resn_sq))
                          + (ee_sq + ec_sq) / v);
      double x_half = x / 2.0;
      // function.cc: tbl_alogch(X) approximates log(cosh(X)).
      double ln_cosh = x_half + std::log((1. + std::exp(-2. * x_half)) / 2.);
      ll = ll_core / 2.0 + ln_cosh;
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

  //! Gradient of the LLGI target for one Miller index w.r.t. sigmaA and
  //! ScatFrac, for use by a sigmaA(resolution)/ScatFrac(resolution)
  //! estimator (see doc/llgi_target_design.md sec. 5) that parameterises
  //! both as B-spline curves and optimises their coefficients directly
  //! against the LLGI target.
  /*! Unlike d_target_one_h_over_fc (whose chain rule through EC = D*fc is
      simple because V does not depend on fc), D = dobs*sigmaa/sqrt(
      scatfrac) enters *both* EC = D*fc *and* V = teps*resn^2*(teps-D^2),
      so this derivative genuinely differs from, and is not obtainable by
      symmetry from, d_target_one_h_over_fc. Derived by hand and verified
      against central finite differences (both centric and acentric,
      max abs error ~1e-9) before being committed here; see also
      tst_llgi.py's exercise_d_target_over_d_sigmaa_scatfrac. Returns
      (d target/d sigmaa, d target/d scatfrac); both already include the
      sign flip matching target_one_h's minimize-me convention.
  */
  inline
  std::pair<double, double>
  d_target_one_h_over_sigmaa_scatfrac(
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
      return std::make_pair(0.0, 0.0);
    }
    double d = dobs * (sigmaa / std::sqrt(scatfrac)) * k;
    double v_e = teps - d * d;
    if (v_e <= 0.0) {
      return std::make_pair(0.0, 0.0);
    }
    double resn_sq = resn * resn;
    double v = teps * resn_sq * v_e;
    double ec = d * fc;
    double ec_sq = ec * ec;
    // dV_by_dD = teps*resn^2 * d(teps - D^2)/dD = -2*D*teps*resn^2
    double dv_by_dd = -2. * d * teps * resn_sq;
    // term1 + term2: derivative of ll_core = -(log(V/(teps^2 resn^2))
    // + (feff^2+EC^2)/V) w.r.t. D.
    double term1 = -(1. / v) * dv_by_dd;
    double term2 = -(2. * ec * fc * v - (feff * feff + ec_sq) * dv_by_dd)
                   / (v * v);
    double d_ll_core_by_dd = term1 + term2;
    double dll_by_dd; // d(ll)/dD, ll in gain form (pre sign-flip)
    if(!centric) {
      double x = 2. * feff * ec / v;
      double bess_term = scitbx::math::bessel::i1_over_i0(x);
      // dX_by_dD = 2*feff*(fc*V - EC*dV_by_dD)/V^2
      double dx_by_dd = 2. * feff * (fc * v - ec * dv_by_dd) / (v * v);
      dll_by_dd = d_ll_core_by_dd + bess_term * dx_by_dd;
    }
    else {
      double x_half = feff * ec / v; // = X/2, matching llgi.h's centric X
      double bess_term = std::tanh(x_half);
      // d(x_half)_by_dD = feff*(fc*V - EC*dV_by_dD)/V^2
      double dxhalf_by_dd = feff * (fc * v - ec * dv_by_dd) / (v * v);
      dll_by_dd = d_ll_core_by_dd / 2.0 + bess_term * dxhalf_by_dd;
    }
    // Sign-flipped to match target_one_h's minimization convention
    // (target_one_h returns -ll).
    double d_target_by_dd = -dll_by_dd;
    // Chain rule through D = dobs*sigmaa*k/sqrt(scatfrac):
    double dd_by_dsigmaa = dobs * k / std::sqrt(scatfrac);
    double dd_by_dscatfrac = -d / (2. * scatfrac);
    double d_target_by_dsigmaa = d_target_by_dd * dd_by_dsigmaa;
    double d_target_by_dscatfrac = d_target_by_dd * dd_by_dscatfrac;
    return std::make_pair(d_target_by_dsigmaa, d_target_by_dscatfrac);
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

  //! Summed LLGI target and per-reflection d(target)/d(sigmaa),
  //! d(target)/d(scatfrac), for a selected set of reflections (in
  //! practice, the R-free/test set -- see doc/llgi_target_design.md sec.
  //! 5.2), for use by a sigmaA(resolution)/ScatFrac(resolution) B-spline
  //! coefficient optimiser. Unlike target_and_gradients (used for the
  //! atomic-parameter refinement target, gradient w.r.t. Fcalc, work-set
  //! normalised), this class:
  //!  - includes every selected reflection in the summed target/gradients
  //!    (no separate work/test split -- the caller already restricts
  //!    `selection` to the desired set, typically r_free_flags itself),
  //!  - returns *per-reflection*, not per-work-set-index, gradient
  //!    arrays (same size and order as the inputs), leaving the chain
  //!    rule through B-spline coefficients (a fixed design matrix, linear
  //!    in the coefficients) to the Python-side estimator,
  //!  - takes sigmaa/scatfrac as already-evaluated per-reflection values
  //!    (i.e. the spline curve evaluated at each reflection's resolution,
  //!    computed Python-side) rather than a resolution parameterisation
  //!    itself, keeping this class's C++ math independent of the spline
  //!    representation.
  class sigmaa_scatfrac_target_and_gradients
  {
    protected:
      double target_;
      af::shared<double> d_target_by_dsigmaa_;
      af::shared<double> d_target_by_dscatfrac_;

    public:
      double target() const { return target_; }
      af::shared<double> const& d_target_by_dsigmaa() const {
        return d_target_by_dsigmaa_;
      }
      af::shared<double> const& d_target_by_dscatfrac() const {
        return d_target_by_dscatfrac_;
      }

      sigmaa_scatfrac_target_and_gradients(
        af::const_ref<double> const& f_eff,
        af::const_ref<bool> const& selection,
        af::const_ref<std::complex<double> > const& f_calc,
        af::const_ref<double> const& dobs,
        af::const_ref<double> const& sigmaa,
        af::const_ref<double> const& scatfrac,
        double scale_factor,
        af::const_ref<double> const& teps,
        af::const_ref<double> const& resn,
        af::const_ref<bool> const& centric_flags)
      :
        target_(0),
        d_target_by_dsigmaa_(f_eff.size(), 0.0),
        d_target_by_dscatfrac_(f_eff.size(), 0.0)
      {
        CCTBX_ASSERT(selection.size() == f_eff.size());
        CCTBX_ASSERT(f_calc.size() == f_eff.size());
        CCTBX_ASSERT(dobs.size() == f_eff.size());
        CCTBX_ASSERT(sigmaa.size() == f_eff.size());
        CCTBX_ASSERT(scatfrac.size() == f_eff.size());
        CCTBX_ASSERT(teps.size() == f_eff.size());
        CCTBX_ASSERT(resn.size() == f_eff.size());
        CCTBX_ASSERT(centric_flags.size() == f_eff.size());
        std::size_t n_selected = 0;
        for(std::size_t i=0;i<f_eff.size();i++) {
          if (!selection[i]) continue;
          n_selected++;
          double feff = f_eff[i];
          double fc = std::abs(f_calc[i]);
          double do_ = dobs[i];
          double sa = sigmaa[i];
          double sf = scatfrac[i];
          double tp = teps[i];
          double rn = resn[i];
          bool c = centric_flags[i];
          target_ += target_one_h(
            feff, fc, do_, sa, sf, scale_factor, tp, rn, c);
          std::pair<double, double> grad = d_target_one_h_over_sigmaa_scatfrac(
            feff, fc, do_, sa, sf, scale_factor, tp, rn, c);
          d_target_by_dsigmaa_[i] = grad.first;
          d_target_by_dscatfrac_[i] = grad.second;
        }
        if (n_selected > 0) {
          double one_over_n = 1. / n_selected;
          target_ *= one_over_n;
          for(std::size_t i=0;i<f_eff.size();i++) {
            d_target_by_dsigmaa_[i] *= one_over_n;
            d_target_by_dscatfrac_[i] *= one_over_n;
          }
        }
      }
  };

}}}} // namespace cctbx::xray::targets::llgi

#endif // GUARD
