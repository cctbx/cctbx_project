#ifndef CCTBX_XRAY_TARGETS_LLGI_E_H
#define CCTBX_XRAY_TARGETS_LLGI_E_H

#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/constants.h>
#include <scitbx/math/bessel.h>
#include <complex>

namespace cctbx { namespace xray { namespace targets { namespace llgi_e {

  //! Log-likelihood-gain-of-intensities target for one miller index,
  //! operating directly on normalised (E-value-scale) amplitudes.
  /*! Companion to cctbx::xray::targets::llgi::target_one_h (llgi.h), which
      is deliberately F-scale (so overall/anisotropic B-factor refinement
      sees Fcalc's resolution falloff directly). This functor exists for a
      different purpose -- see doc/llgi_target_design.md, "An E-Scale LLGI
      Target for Bulk Solvent & SigmaA" -- where normalising out overall
      scale/anisotropy is exactly what is wanted, so the iterated sigmaA/
      bulk-solvent fit doesn't need to wait on or interact with that part
      of the scaling machinery.

      No change of variables is needed here: this is a direct port of
      phasertng's own E-scale functor (phasertng::likelihood::llgi::
      function::target_gradient_hessian, phasertng/codebase/phasertng/
      math/likelihood/llgi/function.cc), specialised to TEPS == 1 (tNCS is
      not yet supported -- see phenix.refinement.llgi_data.
      check_teps_no_tncs, which refuses to load data with TEPS != 1 before
      any of this code ever runs, rather than risk a wrong formula: an
      earlier draft of this derivation proposed V = TEPS*(1-D^2), which
      turned out not to match phasertng's actual V = TEPS - DobsSigaSqr --
      traced to Gfunction::calcRefineTerms, where TEPS (EPSFAC) and the
      amplitude-side tNCS decay factor entangle through a non-trivial,
      HKL-dependent interference sum, not a clean algebraic substitution).

      With TEPS == 1:
        V = 1 - D^2,  D = Dobs*sigmaA,  X = 2*D*Eeff*Emodel/V

      There is no ScatFrac argument here (contrast llgi.h's target_one_h):
      on the E-scale, accounting for how much of the expected scattering
      the model reproduces is the job of how Emodel itself is normalised
      (Emodel = f_model_no_aniso_scale / sqrt(EPS*SigmaP), see the design
      note sec. 4), not of an extra multiplicative term on D.

      eeff     = normalised experimental amplitude, Feff/RESN (RESN is
                 nacelle's own "Root-EpsilonSigmaN" normaliser -- already
                 epsilon- and Wilson-trend-corrected, so no further
                 normalisation of Feff is needed here).
      emodel   = normalised model amplitude, |f_model_no_aniso_scale| /
                 sqrt(EPS * SigmaP) -- see mmtbx.refinement.llgi_sigmaa's
                 Emodel builder.
      dobs     = per-reflection reliability weight (Dobs), from
                 phasertng.nacelle's DOBS column.
      sigmaa   = sigmaA(resolution), a smooth function of resolution only,
                 fitted against this same target on the R-free set (see
                 sigmaa_target_and_gradients below).
      centric  = flag (false for acentric, true for centric).

      Returns a log-likelihood-*gain*, negated so that, like llgi.h's
      target_one_h, it is a quantity for the refinement machinery to
      *minimize*.
  */
  inline
  double
  target_one_h(
    double eeff,
    double emodel,
    double dobs,
    double sigmaa,
    bool centric)
  {
    if (dobs <= 0.0 || sigmaa <= 0.0 || eeff <= 0.0 || emodel <= 0.0) {
      return 0.0;
    }
    double d = dobs * sigmaa;
    double v = 1.0 - d * d;
    if (v <= 0.0) {
      // function.cc: negative_variance() guard; no contribution.
      return 0.0;
    }
    // Symmetrised form. phasertng's function.cc carries the Wilson/null-
    // hypothesis baseline as a separate additive term wll = EOBS^2/teps
    // (teps == 1 here, so wll = eeff^2), added AFTER subtracting
    // eeff^2/V. Those two eeff terms collapse exactly:
    //
    //   -eeff^2/V + eeff^2 = eeff^2*(V - 1)/V = -d^2*eeff^2/V   (V = 1-d^2)
    //
    // so the gain is obtained directly by scaling eeff^2 by d^2 -- exactly
    // the same d^2 factor that turns emodel^2 into ec_sq = (d*emodel)^2.
    // The target is then manifestly SYMMETRIC in eeff and emodel, with no
    // separate baseline term to add and subtract. Verified algebraically
    // and numerically against the original three-term form (agreement to
    // machine precision, both centric and acentric).
    //
    // Note the centric branch works out the same way: halving the core
    // halves -eeff^2/V, which is exactly matched by function.cc halving
    // wll, so the same -d^2*eeff^2/V term serves after halving.
    double d_sq = d * d;
    double ee_sq = d_sq * eeff * eeff;   // symmetric partner of ec_sq
    double ec = d * emodel;              // function.cc: EC = |DobsSigaEcalc|
    double ec_sq = ec * ec;
    double x = 2. * eeff * ec / v;       // function.cc: X
    double ll;
    if(!centric) {
      ll = -(std::log(v) + (ee_sq + ec_sq) / v);
      ll += scitbx::math::bessel::ln_of_i0(x); // function.cc: tbl_alogchI0
    }
    else {
      double ll_core = -(std::log(v) + (ee_sq + ec_sq) / v);
      double x_half = x / 2.0;
      // function.cc: tbl_alogch(X) approximates log(cosh(X)).
      double ln_cosh = x_half + std::log((1. + std::exp(-2. * x_half)) / 2.);
      ll = ll_core / 2.0 + ln_cosh;
    }
    // ll is the log-likelihood-gain; negate to match this file's (and
    // llgi.h's) minimize-me convention.
    return -ll;
  }

  /* \brief Gradient of the E-scale LLGI target for one Miller index
     w.r.t. Emodel (real-valued, since Emodel enters target_one_h only
     through its magnitude here -- the caller, mmtbx.refinement.
     llgi_sigmaa's bulk-solvent chain rule, is responsible for the
     complex-to-real projection through Emodel's phase, exactly as
     llgi.h's d_target_one_h_over_fc does for Fcalc; see that file's
     docstring for the pattern this mirrors). Ported from function.cc's
     dLL_by_dEC, specialised to dEC_by_dEmodel = D (real, positive),
     TEPS == 1.
  */
  inline
  double
  d_target_one_h_over_emodel(
    double eeff,
    double emodel,
    double dobs,
    double sigmaa,
    bool centric)
  {
    if (dobs <= 0.0 || sigmaa <= 0.0 || eeff <= 0.0 || emodel <= 0.0) {
      return 0.0;
    }
    double d = dobs * sigmaa;
    double v = 1.0 - d * d;
    if (v <= 0.0) {
      return 0.0;
    }
    double ec = d * emodel; // function.cc: EC
    double x = 2. * eeff * ec / v; // function.cc: X
    double dll_by_dec; // function.cc: dLL_by_dEC
    if(!centric) {
      double bess_term = scitbx::math::bessel::i1_over_i0(x);
      dll_by_dec = (2. / v) * (eeff * bess_term - ec);
    }
    else {
      double bess_term = std::tanh(x / 2.0);
      dll_by_dec = (1. / v) * (eeff * bess_term - ec);
    }
    double dll_by_demodel = dll_by_dec * d; // dEC_by_dEmodel = d
    // Sign-flipped to match target_one_h's minimization convention.
    return -dll_by_demodel;
  }

  //! Gradient of the E-scale LLGI target for one Miller index w.r.t.
  //! sigmaA, for use by a sigmaA(resolution) estimator (analogous to
  //! llgi.h's d_target_one_h_over_sigmaa_scatfrac, but simpler: there is
  //! no ScatFrac here, and D = dobs*sigmaa enters both EC = D*Emodel and
  //! V = 1-D^2 exactly as in llgi.h, so this shares that derivation's
  //! shape). Derived by hand and verified against central finite
  //! differences (both centric and acentric) before being committed here;
  //! see tst_llgi_e.py's exercise_d_target_over_dsigmaa. Already includes
  //! the sign flip matching target_one_h's minimize-me convention.
  inline
  double
  d_target_one_h_over_sigmaa(
    double eeff,
    double emodel,
    double dobs,
    double sigmaa,
    bool centric)
  {
    if (dobs <= 0.0 || sigmaa <= 0.0 || eeff <= 0.0 || emodel <= 0.0) {
      return 0.0;
    }
    double d = dobs * sigmaa;
    double v = 1.0 - d * d;
    if (v <= 0.0) {
      return 0.0;
    }
    double ec = d * emodel;
    double ec_sq = ec * ec;
    double eeff_sq = eeff * eeff;
    // dV_by_dD = d(1 - D^2)/dD = -2*D
    double dv_by_dd = -2. * d;
    // ll_core = -(log(V) + (Eeff^2+EC^2)/V); dEC_by_dD = Emodel.
    double term1 = -(1. / v) * dv_by_dd;
    double term2 = -(2. * ec * emodel * v - (eeff_sq + ec_sq) * dv_by_dd)
                   / (v * v);
    double d_ll_core_by_dd = term1 + term2;
    double dll_by_dd; // d(ll)/dD, ll in gain form (pre sign-flip)
    if(!centric) {
      double x = 2. * eeff * ec / v;
      double bess_term = scitbx::math::bessel::i1_over_i0(x);
      // dX_by_dD = 2*eeff*(emodel*V - EC*dV_by_dD)/V^2
      double dx_by_dd = 2. * eeff * (emodel * v - ec * dv_by_dd) / (v * v);
      dll_by_dd = d_ll_core_by_dd + bess_term * dx_by_dd;
    }
    else {
      double x_half = eeff * ec / v; // = X/2
      double bess_term = std::tanh(x_half);
      double dxhalf_by_dd = eeff * (emodel * v - ec * dv_by_dd) / (v * v);
      dll_by_dd = d_ll_core_by_dd / 2.0 + bess_term * dxhalf_by_dd;
    }
    double d_target_by_dd = -dll_by_dd; // sign flip, target_one_h convention
    // Chain rule through D = dobs*sigmaa: dd_by_dsigmaa = dobs.
    return d_target_by_dd * dobs;
  }

  //! Summed E-scale LLGI target and per-reflection d(target)/d(sigmaa),
  //! for a selected set of reflections (in practice the R-free/test set --
  //! see doc/llgi_target_design.md's E-scale design note sec. 5), for use
  //! by a sigmaA(resolution) B-spline coefficient optimiser. Mirrors
  //! llgi.h's sigmaa_scatfrac_target_and_gradients: every selected
  //! reflection contributes to the mean (no separate work/test split --
  //! the caller restricts `selection` itself), sigmaa is taken as an
  //! already-evaluated per-reflection value (the spline curve evaluated
  //! at this reflection's resolution, computed Python-side), and the
  //! chain rule through B-spline coefficients is left to that Python-side
  //! estimator.
  class sigmaa_target_and_gradients
  {
    protected:
      double target_;
      af::shared<double> d_target_by_dsigmaa_;

    public:
      double target() const { return target_; }
      af::shared<double> const& d_target_by_dsigmaa() const {
        return d_target_by_dsigmaa_;
      }

      sigmaa_target_and_gradients(
        af::const_ref<double> const& e_eff,
        af::const_ref<bool> const& selection,
        af::const_ref<double> const& e_model,
        af::const_ref<double> const& dobs,
        af::const_ref<double> const& sigmaa,
        af::const_ref<bool> const& centric_flags)
      :
        target_(0),
        d_target_by_dsigmaa_(e_eff.size(), 0.0)
      {
        CCTBX_ASSERT(selection.size() == e_eff.size());
        CCTBX_ASSERT(e_model.size() == e_eff.size());
        CCTBX_ASSERT(dobs.size() == e_eff.size());
        CCTBX_ASSERT(sigmaa.size() == e_eff.size());
        CCTBX_ASSERT(centric_flags.size() == e_eff.size());
        std::size_t n_selected = 0;
        for(std::size_t i=0;i<e_eff.size();i++) {
          if (!selection[i]) continue;
          n_selected++;
          double eeff = e_eff[i];
          double em = e_model[i];
          double do_ = dobs[i];
          double sa = sigmaa[i];
          bool c = centric_flags[i];
          target_ += target_one_h(eeff, em, do_, sa, c);
          d_target_by_dsigmaa_[i] = d_target_one_h_over_sigmaa(
            eeff, em, do_, sa, c);
        }
        if (n_selected > 0) {
          double one_over_n = 1. / n_selected;
          target_ *= one_over_n;
          for(std::size_t i=0;i<e_eff.size();i++) {
            d_target_by_dsigmaa_[i] *= one_over_n;
          }
        }
      }
  };

  //! Summed E-scale LLGI target and per-reflection d(target)/d(Emodel),
  //! for a selected set of reflections (in practice the working set --
  //! all reflections excluding R-free -- see the design note sec. 6), for
  //! use by the bulk-solvent (k_sol, B_sol) chain-rule gradient: the
  //! caller multiplies d_target_by_demodel by d(Emodel)/d(f_model_no_
  //! aniso_scale) (~ 1/sqrt(EPS*SigmaP), Python-side) and then by
  //! d(k_mask*F_mask)/d(k_sol, B_sol) (existing mmtbx.bulk_solvent
  //! derivative code) to get d(target)/d(k_sol, B_sol). sigmaA here is
  //! fixed (this stage's whole point is to hold sigmaA(d) constant while
  //! refitting bulk solvent -- see design note sec. 7).
  class emodel_target_and_gradients
  {
    protected:
      double target_;
      af::shared<double> d_target_by_demodel_;

    public:
      double target() const { return target_; }
      af::shared<double> const& d_target_by_demodel() const {
        return d_target_by_demodel_;
      }

      emodel_target_and_gradients(
        af::const_ref<double> const& e_eff,
        af::const_ref<bool> const& selection,
        af::const_ref<double> const& e_model,
        af::const_ref<double> const& dobs,
        af::const_ref<double> const& sigmaa,
        af::const_ref<bool> const& centric_flags)
      :
        target_(0),
        d_target_by_demodel_(e_eff.size(), 0.0)
      {
        CCTBX_ASSERT(selection.size() == e_eff.size());
        CCTBX_ASSERT(e_model.size() == e_eff.size());
        CCTBX_ASSERT(dobs.size() == e_eff.size());
        CCTBX_ASSERT(sigmaa.size() == e_eff.size());
        CCTBX_ASSERT(centric_flags.size() == e_eff.size());
        std::size_t n_selected = 0;
        for(std::size_t i=0;i<e_eff.size();i++) {
          if (!selection[i]) continue;
          n_selected++;
          double eeff = e_eff[i];
          double em = e_model[i];
          double do_ = dobs[i];
          double sa = sigmaa[i];
          bool c = centric_flags[i];
          target_ += target_one_h(eeff, em, do_, sa, c);
          d_target_by_demodel_[i] = d_target_one_h_over_emodel(
            eeff, em, do_, sa, c);
        }
        if (n_selected > 0) {
          double one_over_n = 1. / n_selected;
          target_ *= one_over_n;
          for(std::size_t i=0;i<e_eff.size();i++) {
            d_target_by_demodel_[i] *= one_over_n;
          }
        }
      }
  };

}}}} // namespace cctbx::xray::targets::llgi_e

#endif // CCTBX_XRAY_TARGETS_LLGI_E_H
