/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     Oct 2002: Modified copy of phenix/fast_translation/combinations.h (rwgk)
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_TRANSLATION_SEARCH_FAST_NV1995_COMBINATIONS_H
#define CCTBX_TRANSLATION_SEARCH_FAST_NV1995_COMBINATIONS_H

#include <cctbx/translation_search/fast_nv1995/intermediates.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace cctbx { namespace translation_search { namespace fast_nv1995_detail {

  // Navaza & Vernoslova (1995), p. 447, Eq. (13)
  template <typename FloatTypeSummation,
            typename FloatTypeTarget>
  void
  combination_eq13(
    intermediates<FloatTypeSummation> const& interm,
    af::const_ref<FloatTypeTarget, af::c_grid<3> > const&
      result_eq14,
    af::ref<FloatTypeTarget, af::c_grid<3> > const&
      target_accu) // on input result of equation 15
  {
    CCTBX_ASSERT(        result_eq14.accessor()
                 .all_eq(target_accu.accessor()));
    typedef FloatTypeTarget f_t;
    for (std::size_t i=0;i<target_accu.size();i++) {
      f_t r14p = result_eq14[i] + interm.sum_mh_f_sq;
      f_t r15p = target_accu[i] + interm.sum_mh_f_sq_f_sq;
      f_t d = r15p - r14p * r14p / interm.sum_mh;
      if (d > f_t(0)) {
        target_accu[i] = std::sqrt(d) * std::sqrt(interm.sum_mh_d_i_obs_sq);
      }
      else {
        target_accu[i] = f_t(0);
      }
    }
  }

  // Navaza & Vernoslova (1995), p. 447, Eq. (13)
  template <typename FloatTypeSummation,
            typename FloatTypeTarget>
  void
  combination_eq12(
    intermediates<FloatTypeSummation> const& interm,
    af::const_ref<FloatTypeTarget, af::c_grid<3> > const&
      result_eq14_with_i_obs,
    af::ref<FloatTypeTarget, af::c_grid<3> > const&
      target_accu, // on input result of equation 13
    FloatTypeTarget const& big_correlation=1.e6)
  {
    CCTBX_ASSERT(        result_eq14_with_i_obs.accessor()
                 .all_eq(target_accu.accessor()));
    typedef FloatTypeTarget f_t;
    for (std::size_t i=0;i<target_accu.size();i++) {
      f_t r14ip = result_eq14_with_i_obs[i] + interm.sum_mh_f_sq_d_i_obs;
      if (scitbx::fn::absolute(r14ip / big_correlation) < target_accu[i]) {
        target_accu[i] = r14ip / target_accu[i];
      }
      else {
        target_accu[i] = f_t(0);
      }
    }
  }

}}} // namespace cctbx::translation_search::fast_nv1995_detail

#endif // CCTBX_TRANSLATION_SEARCH_FAST_NV1995_COMBINATIONS_H
