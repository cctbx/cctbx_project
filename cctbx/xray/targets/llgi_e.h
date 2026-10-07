#ifndef CCTBX_XRAY_TARGETS_LLGI_E_H
#define CCTBX_XRAY_TARGETS_LLGI_E_H

#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/shared.h>
#include <cctbx/xray/targets/llgi.h>
#include <cctbx/xray/targets/llgi_exact.h>

namespace cctbx { namespace xray { namespace targets { namespace llgi_e {

  //! Log-likelihood-gain-of-intensities target for one miller index,
  //! on normalised (E-value-scale) amplitudes.
  /*! The E-scale counterpart of llgi::target_one_h (llgi.h), used for the
      sigmaA(resolution) fit, where normalising out the overall scale and
      anisotropy is what is wanted. With TEPS == 1 it is llgi.h's target
      with ScatFrac = RESN = TEPS = k = 1:
        V = 1 - D^2,  D = Dobs*sigmaA,  X = 2*D*Eeff*Emodel/V
      There is no ScatFrac: on the E scale, the fraction of the expected
      scattering the model accounts for is absorbed by normalising Emodel
      (mmtbx.refinement.llgi_e_bulk_solvent.build_e_model).

      eeff     = Feff/RESN.
      emodel   = |f_model_no_aniso_scale|/sqrt(EPS*SigmaP).
      dobs     = Dobs (nacelle DOBS column).
      sigmaa   = sigmaA(resolution).
      centric  = flag (false for acentric, true for centric).

      Returns the negated log-likelihood gain (minimize-me convention).
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
    return llgi::target_one_h(
      eeff, emodel, dobs, sigmaa, 1., 1., 1., 1., centric);
  }

  //! Gradient of target_one_h w.r.t. sigmaA (minimize-me convention).
  inline
  double
  d_target_one_h_over_sigmaa(
    double eeff,
    double emodel,
    double dobs,
    double sigmaa,
    bool centric)
  {
    return llgi::d_target_one_h_over_sigmaa_scatfrac(
      eeff, emodel, dobs, sigmaa, 1., 1., 1., 1., centric).first;
  }

  //! Mean E-scale LLGI target and per-reflection d(target)/d(sigmaa) over
  //! a selected set of reflections (the sigmaA fit uses the R-free set),
  //! as llgi::sigmaa_scatfrac_target_and_gradients does on the F scale:
  //! sigmaa is an already-evaluated per-reflection value, and the chain
  //! rule through the curve's parameters is left to the Python-side fit.
  class sigmaa_target_and_gradients
  {
    protected:
      double target_;
      af::shared<double> d_target_by_dsigmaa_;
      std::size_t n_exact_;

    public:
      double target() const { return target_; }
      af::shared<double> const& d_target_by_dsigmaa() const {
        return d_target_by_dsigmaa_;
      }
      //! Number of selected reflections evaluated with the exact LLGI.
      std::size_t n_exact() const { return n_exact_; }

      //! hybrid (optional): use the exact LLGI where it says so (see
      //! llgi_exact::hybrid); its arrays are indexed like e_eff.
      sigmaa_target_and_gradients(
        af::const_ref<double> const& e_eff,
        af::const_ref<bool> const& selection,
        af::const_ref<double> const& e_model,
        af::const_ref<double> const& dobs,
        af::const_ref<double> const& sigmaa,
        af::const_ref<bool> const& centric_flags,
        llgi_exact::hybrid const* hybrid = 0)
      :
        target_(0),
        d_target_by_dsigmaa_(e_eff.size(), 0.0),
        n_exact_(0)
      {
        CCTBX_ASSERT(hybrid == 0 || hybrid->size() == e_eff.size());
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
          if (hybrid != 0 && hybrid->use_exact(i, sa)) {
            llgi_exact::result r = hybrid->evaluate_at(i, em, sa, c);
            target_ -= r.ll;
            d_target_by_dsigmaa_[i] = -r.d_ll_d_a;
            n_exact_++;
            continue;
          }
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

}}}} // namespace cctbx::xray::targets::llgi_e

#endif // CCTBX_XRAY_TARGETS_LLGI_E_H
