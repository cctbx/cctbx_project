#ifndef CCTBX_XRAY_TARGETS_LLGI_H
#define CCTBX_XRAY_TARGETS_LLGI_H

#include <cctbx/xray/targets.h>
#include <scitbx/constants.h>
#include <scitbx/math/bessel.h>
#include <utility>
#include <cctbx/xray/targets/llgi_exact.h>

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
        Ecalc = k*Fcalc/sqrt(TEPS*RESN^2*ScatFrac)     -> fc   (unnormalized,
                                                                 |Fcalc|)
        D = Dobs*sigmaA (mean D*Ecalc)                 -> d = D*k/sqrt(
                                                                 ScatFrac)
        V_E = TEPS - D^2                               -> V_F = TEPS*RESN^2
                                                                 *V_E
      ScatFrac is the model's own normalisation: the mean of
      k^2|Fcalc|^2/(TEPS*RESN^2) as a function of resolution (the fraction
      of the expected intensity the model accounts for), so that Ecalc has
      unit mean square. It enters the mean only; the variance depends on
      sigmaA alone, as in phasertng (V = TEPS - DobsSigaSqr).

      feff     = effective amplitude (Feff), NOT normalized to the E-value
                 scale. Supplied by phasertng.nacelle as the FEFF column.
      fc       = |Fcalc|, NOT normalized (ordinary calculated structure
                 factor amplitude on the same scale as feff).
      dobs     = per-reflection reliability weight (Dobs), varies from
                 reflection to reflection. Supplied by phasertng.nacelle as
                 the DOBS column.
      sigmaa   = sigmaA(resolution), a smooth function of resolution only.
      scatfrac = ScatFrac(resolution), a smooth function of resolution
                 only (see above).
      k        = overall scale coefficient (values <= 0 are taken as 1).
      teps     = per-reflection tNCS-adjusted normalization factor (nacelle
                 TEPS column); equals 1 in the ordinary case (no tNCS), and
                 must be all ones when the tncs_u/tncs_t arrays are given
                 (the tNCS modulation then comes from tncs_t).
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
    // D = Dobs*sigmaA (function.cc's DobsSiga) sets the variance; the mean
    // of fc is d*fc, d = D*k/sqrt(ScatFrac), from Ecalc's normalisation.
    double big_d = dobs * sigmaa;
    double d = big_d * k / std::sqrt(scatfrac);
    double resn_sq = resn * resn;
    // v_e is the E-scale variance, function.cc's V = teps - DobsSigaSqr;
    // v is v_e rescaled to the F-scale, v = teps*resn^2*v_e.
    double v_e = teps - big_d * big_d;
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
    //   -feff^2/V + feff^2/(teps^2*resn^2) = -(D^2/teps)*feff^2/V
    //
    // (using V = teps*resn^2*(teps - D^2)), so the gain follows directly
    // by scaling feff^2 -- no separate baseline term to add and subtract.
    //
    // The multiplier is D^2/TEPS, not D^2: the two coincide only at
    // TEPS == 1 (tst_llgi.py checks the TEPS != 1 case against the
    // three-term form).
    double ee_sq = (big_d * big_d / teps) * feff_sq;
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
     d = dobs*sigmaa*k/sqrt(scatfrac), so EC = d*fc reduces the chain rule to
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
    double fc = std::abs(fc_complex);
    if(k <= 0.0) k = 1.0;
    std::complex<double> d_target_over_fc(0, 0);
    // Same "no contribution" conditions as target_one_h.
    if (dobs <= 0.0 || sigmaa <= 0.0 || scatfrac <= 0.0
        || feff <= 0.0 || fc <= 0.0) {
      return d_target_over_fc;
    }
    double big_d = dobs * sigmaa;
    double d = big_d * k / std::sqrt(scatfrac);
    double v_e = teps - big_d * big_d;
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
  //! ScatFrac, for the sigmaA/ScatFrac fits.
  /*! sigmaA enters the variance through D = dobs*sigmaa (V = teps*resn^2*
      (teps - D^2), and the data term (D^2/teps)*feff^2) and the mean
      through d = D*k/sqrt(scatfrac) (EC = d*fc); ScatFrac enters the mean
      only. Returns (d target/d sigmaa, d target/d scatfrac), both with the
      sign flip of target_one_h's minimize-me convention (checked by
      finite differences in tst_llgi.py's
      exercise_sigmaa_scatfrac_finite_difference_gradients).
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
    double big_d = dobs * sigmaa;
    double d = big_d * k / std::sqrt(scatfrac);
    double v_e = teps - big_d * big_d;
    if (v_e <= 0.0) {
      return std::make_pair(0.0, 0.0);
    }
    double resn_sq = resn * resn;
    double v = teps * resn_sq * v_e;
    double feff_sq = feff * feff;
    double ec = d * fc;
    double ee_sq = (big_d * big_d / teps) * feff_sq;
    double ec_sq = ec * ec;
    double x = 2. * feff * ec / v;
    // ll_core = -(log(v/(teps^2 resn^2)) + (ee_sq + ec_sq)/v), and the
    // Bessel (acentric) or log cosh (centric, with ll_core/2) term in x.
    // Partial derivatives w.r.t. D (through v and ee_sq) and d (ec):
    double dv_by_dbig_d = -2. * big_d * teps * resn_sq;
    double dee_by_dbig_d = 2. * big_d * feff_sq / teps;
    double dcore_by_dbig_d = -dv_by_dbig_d / v
      - (dee_by_dbig_d * v - (ee_sq + ec_sq) * dv_by_dbig_d) / (v * v);
    double dx_by_dbig_d = -x * dv_by_dbig_d / v;
    double dcore_by_dd = -2. * d * fc * fc / v;
    double dx_by_dd = 2. * feff * fc / v;
    double dll_by_dbig_d, dll_by_dd; // ll in gain form (pre sign-flip)
    if(!centric) {
      double bess_term = scitbx::math::bessel::i1_over_i0(x);
      dll_by_dbig_d = dcore_by_dbig_d + bess_term * dx_by_dbig_d;
      dll_by_dd = dcore_by_dd + bess_term * dx_by_dd;
    }
    else {
      double tanh_term = std::tanh(x / 2.0);
      dll_by_dbig_d = dcore_by_dbig_d / 2.0 + tanh_term * dx_by_dbig_d / 2.0;
      dll_by_dd = dcore_by_dd / 2.0 + tanh_term * dx_by_dd / 2.0;
    }
    // D = dobs*sigmaa, d = dobs*sigmaa*k/sqrt(scatfrac)
    double dll_by_dsigmaa = dll_by_dbig_d * dobs
                          + dll_by_dd * dobs * k / std::sqrt(scatfrac);
    double dll_by_dscatfrac = dll_by_dd * (-d / (2. * scatfrac));
    // Sign-flipped to match target_one_h's minimization convention.
    return std::make_pair(-dll_by_dsigmaa, -dll_by_dscatfrac);
  }

  //! Ecalc/|Fcalc| for the exact LLGI and the tNCS frame (TEPS == 1):
  //! Ecalc = k*|Fcalc|/(RESN*sqrt(ScatFrac)), the unit-normalised model
  //! amplitude, whose coefficient in the likelihood is sigmaA itself.
  //! 0 where ScatFrac or RESN is not positive.
  inline
  double
  ecalc_per_fc(double scatfrac, double k, double resn)
  {
    if (scatfrac <= 0.0 || resn <= 0.0) return 0.0;
    if (k <= 0.0) k = 1.0;
    return k / (resn * std::sqrt(scatfrac));
  }

  //! Whether the hybrid LLGI uses the exact likelihood for reflection i.
  inline
  bool
  use_exact(
    llgi_exact::hybrid const* hybrid, std::size_t i, double sigmaa,
    double teps, double resn)
  {
    return hybrid != 0 && std::abs(teps - 1.0) < 1e-6 && resn > 0
        && hybrid->use_exact(i, sigmaa);
  }

  //! Translational NCS (tNCS) on the E scale normalised by RESN (without
  //! the tNCS modulation): the model term has mean a*ec and the error
  //! variance is (1 - a^2)*u, where u is the tNCS modulation of the error
  //! variance; the data are normalised by t, the tNCS modulation of the
  //! expected intensity, which must exceed (1 - a^2)*u. In the
  //! t-normalised frame this is the ordinary LLGI (variance 1 - a'^2 about
  //! a'*ec') with
  //!   a' = sqrt(1 - (1 - a^2)*u/t),   ec' = s*ec,   s = a/(sqrt(t)*a').
  //! With u = t = 1, a' = a and s = 1. valid is false (no contribution)
  //! unless 0 < a < 1, u > 0 and t > (1 - a^2)*u.
  struct tncs_frame
  {
    bool valid;
    double a_p;      // a'
    double s;        // ec'/ec
    double da_p_da;  // da'/da
    double ds_da;    // ds/da

    tncs_frame(double a, double u, double t)
    :
      valid(false), a_p(0), s(0), da_p_da(0), ds_da(0)
    {
      if (!(a > 0 && a < 1 && u > 0 && t > 0)) return;
      double a_p_sq = 1. - (1. - a * a) * u / t;
      if (!(a_p_sq > 0)) return;
      a_p = std::sqrt(a_p_sq);
      double sqrt_t = std::sqrt(t);
      s = a / (sqrt_t * a_p);
      da_p_da = a * u / (t * a_p);
      ds_da = (1. - a * da_p_da / a_p) / (sqrt_t * a_p);
      valid = true;
    }
  };

  //! tncs_frame for arrays of a (sigmaA), u and t: a' and s per reflection,
  //! 0 where the frame is invalid (no contribution), and the valid flags.
  struct tncs_frame_many
  {
    af::shared<double> a_p, s;
    af::shared<bool> valid;

    tncs_frame_many(
      af::const_ref<double> const& a,
      af::const_ref<double> const& u,
      af::const_ref<double> const& t)
    {
      CCTBX_ASSERT(u.size() == a.size() && t.size() == a.size());
      a_p.reserve(a.size()); s.reserve(a.size()); valid.reserve(a.size());
      for (std::size_t i = 0; i < a.size(); i++) {
        tncs_frame f(a[i], u[i], t[i]);
        a_p.push_back(f.a_p); s.push_back(f.s); valid.push_back(f.valid);
      }
    }
  };

  //! Target (minimize-me) for one reflection in the t-normalised frame of
  //! tncs_frame, and its derivatives w.r.t. ec' and a': the exact LLGI
  //! where the hybrid says so, otherwise the Rice form (target_one_h with
  //! RESN = TEPS = ScatFrac = k = 1). eeff_t = Feff/(RESN*sqrt(t)).
  struct tncs_evaluation
  {
    double target, d_target_d_ec, d_target_d_a;
    bool exact;
  };

  inline
  tncs_evaluation
  evaluate_tncs(
    llgi_exact::hybrid const* hybrid,
    std::size_t i,
    double eeff_t,
    double ec_p,
    double dobs,
    double a_p,
    bool centric,
    bool need_d_a)
  {
    tncs_evaluation e;
    e.d_target_d_a = 0;
    if (hybrid != 0 && hybrid->use_exact(i, a_p)) {
      llgi_exact::result r = hybrid->evaluate_at(i, ec_p, a_p, centric);
      e.target = -r.ll;
      e.d_target_d_ec = -r.d_ll_d_ec;
      e.d_target_d_a = -r.d_ll_d_a;
      e.exact = true;
      return e;
    }
    e.target = target_one_h(eeff_t, ec_p, dobs, a_p, 1., 1., 1., 1., centric);
    e.d_target_d_ec = d_target_one_h_over_fc(eeff_t,
      std::complex<double>(ec_p, 0.), dobs, a_p, 1., 1., 1., 1., centric).real();
    if (need_d_a) {
      e.d_target_d_a = d_target_one_h_over_sigmaa_scatfrac(
        eeff_t, ec_p, dobs, a_p, 1., 1., 1., 1., centric).first;
    }
    e.exact = false;
    return e;
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

      With tNCS (tncs_u, tncs_t non-empty; teps must then be 1): see
      tncs_frame. f_eff and the hybrid's observations must be normalised
      by tncs_t (f_eff = Eeff*RESN*sqrt(t), E_obs^2 = I/(RESN^2*t)).
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
        bool compute_gradients,
        llgi_exact::hybrid const* hybrid = 0)
      :
        common_results(f_eff.size()),
        n_exact_(0)
      {
        compute(f_eff, r_free_flags, f_calc, dobs, sigmaa, scatfrac,
          scale_factor, teps, resn, centric_flags, compute_gradients, hybrid,
          af::const_ref<double>(0, 0), af::const_ref<double>(0, 0));
      }

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
        bool compute_gradients,
        llgi_exact::hybrid const* hybrid,
        af::const_ref<double> const& tncs_u,
        af::const_ref<double> const& tncs_t)
      :
        common_results(f_eff.size()),
        n_exact_(0)
      {
        compute(f_eff, r_free_flags, f_calc, dobs, sigmaa, scatfrac,
          scale_factor, teps, resn, centric_flags, compute_gradients, hybrid,
          tncs_u, tncs_t);
      }

      //! Number of working-set reflections evaluated with the exact LLGI.
      std::size_t n_exact() const { return n_exact_; }

    protected:
      std::size_t n_exact_;

      void
      compute(
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
        bool compute_gradients,
        llgi_exact::hybrid const* hybrid,
        af::const_ref<double> const& tncs_u,
        af::const_ref<double> const& tncs_t)
      {
        CCTBX_ASSERT(tncs_t.size() == tncs_u.size());
        bool tncs = tncs_u.size() != 0;
        CCTBX_ASSERT(!tncs || (tncs_u.size() == f_eff.size()
                            && tncs_t.size() == f_eff.size()));
        CCTBX_ASSERT(hybrid == 0 || hybrid->size() == f_eff.size());
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
          bool work = rffs.is_work_refl(i);
          double ec_per_fc = ecalc_per_fc(sf, scale_factor, rn);
          double t = 0;
          std::complex<double> g(0, 0);
          bool exact = false;
          if (tncs) {
            CCTBX_ASSERT(std::abs(tp - 1.0) < 1e-6);
            tncs_frame fr(sa, tncs_u[i], tncs_t[i]);
            if (fr.valid && ec_per_fc > 0) {
              tncs_evaluation e = evaluate_tncs(hybrid, i,
                feff / (rn * std::sqrt(tncs_t[i])), fr.s * fc * ec_per_fc,
                do_, fr.a_p, c, false);
              t = e.target;
              exact = e.exact;
              if (work && compute_gradients && fc > 0) {
                // ec' = s*Ecalc = s*fc*ec_per_fc
                g = (e.d_target_d_ec * fr.s * ec_per_fc)
                  * (std::conj(f_calc[i]) / fc);
              }
            }
          }
          else {
            exact = ec_per_fc > 0 && use_exact(hybrid, i, sa, tp, rn);
            if (exact) {
              llgi_exact::result r = hybrid->evaluate_at(
                i, fc * ec_per_fc, sa, c);
              t = -r.ll;
              if (work && compute_gradients && fc > 0) {
                // d(-LLGI)/dfc = -(dLLGI/dEcalc)*ec_per_fc along the phase
                g = (-r.d_ll_d_ec * ec_per_fc) * (std::conj(f_calc[i]) / fc);
              }
            }
            else {
              t = target_one_h(feff, fc, do_, sa, sf, scale_factor, tp, rn, c);
              if (work && compute_gradients) {
                g = d_target_one_h_over_fc(
                  feff, f_calc[i], do_, sa, sf, scale_factor, tp, rn, c);
              }
            }
          }
          if (exact && work) n_exact_++;
          target_per_reflection_[i] = t;
          if (work) {
            target_work += t;
            if (compute_gradients) {
              gradients_work_.push_back(std::conj(g) * one_over_n_work);
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

  //! Mean LLGI target and per-reflection d(target)/d(sigmaa),
  //! d(target)/d(scatfrac) over a selected set of reflections (the
  //! ScatFrac fit uses the working set), for the sigmaA/ScatFrac fits.
  //! Unlike target_and_gradients (the atomic-parameter target, gradient
  //! w.r.t. Fcalc, work-set normalised), this class:
  //!  - averages over the selected reflections only (the caller chooses
  //!    the set),
  //!  - returns per-reflection gradient arrays (same size and order as
  //!    the inputs), leaving the chain rule through the curve's parameters
  //!    to the Python-side fit,
  //!  - takes sigmaa/scatfrac as already-evaluated per-reflection values.
  class sigmaa_scatfrac_target_and_gradients
  {
    protected:
      double target_;
      af::shared<double> d_target_by_dsigmaa_;
      af::shared<double> d_target_by_dscatfrac_;
      std::size_t n_exact_;

    public:
      double target() const { return target_; }
      //! Number of selected reflections evaluated with the exact LLGI.
      std::size_t n_exact() const { return n_exact_; }
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
        af::const_ref<bool> const& centric_flags,
        llgi_exact::hybrid const* hybrid = 0)
      :
        target_(0),
        d_target_by_dsigmaa_(f_eff.size(), 0.0),
        d_target_by_dscatfrac_(f_eff.size(), 0.0),
        n_exact_(0)
      {
        compute(f_eff, selection, f_calc, dobs, sigmaa, scatfrac,
          scale_factor, teps, resn, centric_flags, hybrid,
          af::const_ref<double>(0, 0), af::const_ref<double>(0, 0));
      }

      //! With tNCS: see target_and_gradients.
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
        af::const_ref<bool> const& centric_flags,
        llgi_exact::hybrid const* hybrid,
        af::const_ref<double> const& tncs_u,
        af::const_ref<double> const& tncs_t)
      :
        target_(0),
        d_target_by_dsigmaa_(f_eff.size(), 0.0),
        d_target_by_dscatfrac_(f_eff.size(), 0.0),
        n_exact_(0)
      {
        compute(f_eff, selection, f_calc, dobs, sigmaa, scatfrac,
          scale_factor, teps, resn, centric_flags, hybrid, tncs_u, tncs_t);
      }

    protected:
      void
      compute(
        af::const_ref<double> const& f_eff,
        af::const_ref<bool> const& selection,
        af::const_ref<std::complex<double> > const& f_calc,
        af::const_ref<double> const& dobs,
        af::const_ref<double> const& sigmaa,
        af::const_ref<double> const& scatfrac,
        double scale_factor,
        af::const_ref<double> const& teps,
        af::const_ref<double> const& resn,
        af::const_ref<bool> const& centric_flags,
        llgi_exact::hybrid const* hybrid,
        af::const_ref<double> const& tncs_u,
        af::const_ref<double> const& tncs_t)
      {
        CCTBX_ASSERT(tncs_t.size() == tncs_u.size());
        bool tncs = tncs_u.size() != 0;
        CCTBX_ASSERT(!tncs || (tncs_u.size() == f_eff.size()
                            && tncs_t.size() == f_eff.size()));
        CCTBX_ASSERT(hybrid == 0 || hybrid->size() == f_eff.size());
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
          double ec_per_fc = ecalc_per_fc(sf, scale_factor, rn);
          double ec = fc * ec_per_fc; // Ecalc = k*fc/(resn*sqrt(scatfrac))
          if (tncs) {
            CCTBX_ASSERT(std::abs(tp - 1.0) < 1e-6);
            tncs_frame fr(sa, tncs_u[i], tncs_t[i]);
            if (!fr.valid || !(ec_per_fc > 0)) continue;
            tncs_evaluation e = evaluate_tncs(hybrid, i,
              feff / (rn * std::sqrt(tncs_t[i])), fr.s * ec, do_, fr.a_p, c,
              true);
            target_ += e.target;
            // sigmaA through a' and ec' = s*Ecalc; ScatFrac through Ecalc
            d_target_by_dsigmaa_[i] = e.d_target_d_a * fr.da_p_da
                                    + e.d_target_d_ec * fr.ds_da * ec;
            d_target_by_dscatfrac_[i] = e.d_target_d_ec * fr.s
                                      * (-ec / (2. * sf));
            if (e.exact) n_exact_++;
            continue;
          }
          if (ec_per_fc > 0 && use_exact(hybrid, i, sa, tp, rn)) {
            llgi_exact::result r = hybrid->evaluate_at(i, ec, sa, c);
            target_ -= r.ll;
            d_target_by_dsigmaa_[i] = -r.d_ll_d_a;
            // dEcalc/dScatFrac = -Ecalc/(2 ScatFrac)
            d_target_by_dscatfrac_[i] = r.d_ll_d_ec * ec / (2. * sf);
            n_exact_++;
            continue;
          }
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
