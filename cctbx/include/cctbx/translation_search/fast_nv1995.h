/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     Oct 2002: Modified copy of phenix/fast_translation/front_end.h (rwgk)
     Mar 2002: using cctbx/array_family (rwgk)
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_TRANSLATION_SEARCH_FAST_NV1995_H
#define CCTBX_TRANSLATION_SEARCH_FAST_NV1995_H

#include <cctbx/translation_search/fast_nv1995/combinations.h>
#include <cctbx/translation_search/fast_nv1995/summations.h>
#include <cctbx/translation_search/fast_nv1995/fft.h>
#include <cctbx/translation_search/map_gridding.h>
#include <cctbx/miller/index_span.h>
#include <scitbx/fftpack/real_to_complex_3d.h> // XXX n_complex_from_n_real
#include <scitbx/array_family/versa.h>

// Navaza, J. & Vernoslova, E. (1995). Acta Cryst. A51, 445-449.

namespace cctbx { namespace translation_search {

  template <typename FloatTypeSummation = double,
            typename FloatTypeTarget = float>
  class fast_nv1995
  {
    public:
      template <typename FloatTypeFobs,
                typename FloatTypeFpart,
                typename FloatTypeP1Fcalc>
      fast_nv1995(
        map_gridding<> const& gridding,
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<miller::index<> > const& miller_indices_f_obs,
        af::const_ref<FloatTypeFobs> const& f_obs,
        af::const_ref<std::complex<FloatTypeFpart> > const& f_part,
        af::const_ref<miller::index<> > const& miller_indices_p1_f_calc,
        af::const_ref<std::complex<FloatTypeP1Fcalc> > const& p1_f_calc)
      {
        // FUTURE: move out of class body
        using namespace fast_nv1995_detail;
        typedef FloatTypeSummation sum_f_t;
        typedef FloatTypeTarget target_f_t;
        typedef fast_nv1995_detail::summation_accumulator<
                  sum_f_t> sum_accu_type;
        typedef fast_nv1995_detail::shrinking_complex_to_real_fft<
                  sum_f_t, target_f_t> shrinking_fft;

        intermediates<sum_f_t> interm(
          space_group, anomalous_flag, miller_indices_f_obs, f_obs, f_part);

        af::int3 range = miller::index_span(
          miller_indices_p1_f_calc).abs_range();
        f_calc_map<sum_f_t> fc_map(anomalous_flag, range);
        fc_map.import(miller_indices_p1_f_calc, p1_f_calc);

        af::c_grid<3> dim_target(
          gridding.target());
        af::c_grid<3> dim_eq14(
          scitbx::fftpack::n_complex_from_n_real(gridding.target()));
        af::c_grid<3> dim_eq15(
          scitbx::fftpack::n_complex_from_n_real(gridding.target()));

        // separate scope for allocation of large memory block
        {
          af::versa<std::complex<sum_f_t>, af::c_grid<3> > sum_eq15(dim_eq15);
          sum_accu_type accu(
            sum_eq15.begin(), gridding.target(), miller::index<>(dim_eq15));
          summation_eq15(space_group, miller_indices_f_obs,
            f_part, fc_map, accu);
          af::ref<target_f_t, af::c_grid<3> > res_eq15 =
            shrinking_fft::run(sum_eq15.ref(), gridding.target(), dim_target);
          target_map_.as_base_array().assign(res_eq15);
          CCTBX_ASSERT(   target_map_.as_base_array().size()
                       == dim_target.size_1d());
          target_map_.resize(dim_target);
        }

        af::versa<std::complex<sum_f_t>, af::c_grid<3> > sum_eq14(dim_eq14);
        sum_accu_type accu(
          sum_eq14.begin(), gridding.target(), miller::index<>(dim_eq14));
        summation_eq14(space_group, miller_indices_f_obs,
          af::const_ref<sum_f_t>(0,0), f_part, fc_map, accu);
        af::ref<target_f_t, af::c_grid<3> > res_eq14 =
          shrinking_fft::run(sum_eq14.ref(), gridding.target(), dim_target);
        combination_eq13(interm, res_eq14, target_map_.ref());

        sum_eq14.fill(0);
        summation_eq14(space_group, miller_indices_f_obs,
          interm.d_i_obs.const_ref(), f_part, fc_map, accu);
        res_eq14 =
          shrinking_fft::run(sum_eq14.ref(), gridding.target(), dim_target);
        combination_eq12(interm, res_eq14, target_map_.ref());
      }

      af::versa<FloatTypeTarget, af::c_grid<3> >
      target_map() const { return target_map_; }

    protected:
      af::versa<FloatTypeTarget, af::c_grid<3> > target_map_;
  };

}} // namespace cctbx::translation_search

#endif // CCTBX_TRANSLATION_SEARCH_FAST_NV1995_H
