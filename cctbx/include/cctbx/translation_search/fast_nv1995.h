#ifndef CCTBX_TRANSLATION_SEARCH_FAST_NV1995_H
#define CCTBX_TRANSLATION_SEARCH_FAST_NV1995_H

#include <cctbx/translation_search/fast_nv1995/combinations.h>
#include <cctbx/translation_search/fast_nv1995/summations.h>
#include <cctbx/maptbx/copy.h>
#include <cctbx/miller/index_span.h>
#include <scitbx/fftpack/real_to_complex_3d.h>

// Navaza, J. & Vernoslova, E. (1995). Acta Cryst. A51, 445-449.

namespace cctbx { namespace translation_search {

  template <typename FloatType = double>
  class fast_nv1995
  {
    public:
      template <typename FloatTypeFobs,
                typename FloatTypeFpart,
                typename FloatTypeP1Fcalc>
      fast_nv1995(
        af::int3 const& gridding,
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<miller::index<> > const& miller_indices_f_obs,
        af::const_ref<FloatTypeFobs> const& f_obs,
        af::const_ref<std::complex<FloatTypeFpart> > const& f_part,
        af::const_ref<miller::index<> > const& miller_indices_p1_f_calc,
        af::const_ref<std::complex<FloatTypeP1Fcalc> > const& p1_f_calc)
      :
        target_map_(af::c_grid<3>(gridding))
      {
        // FUTURE: move out of class body
        using namespace fast_nv1995_detail;
        typedef FloatType f_t;

        intermediates<f_t> interm(
          space_group, anomalous_flag, miller_indices_f_obs, f_obs, f_part);

        af::int3 range = miller::index_span(
          miller_indices_p1_f_calc).abs_range();
        f_calc_map<f_t> fc_map(anomalous_flag, range);
        fc_map.import(miller_indices_p1_f_calc, p1_f_calc);

        scitbx::fftpack::real_to_complex_3d<f_t> rfft(gridding);

        af::versa<std::complex<f_t>, af::c_grid<3> >
          accu_mem_complex(
            af::c_grid<3>(rfft.n_complex()));
        af::const_ref<f_t, af::c_grid_padded<3> >
          accu_mem_real_const_ref(
            reinterpret_cast<f_t*>(accu_mem_complex.begin()),
            af::c_grid_padded<3>(rfft.m_real(), rfft.n_real()));

        summation_accumulator<f_t>
          accu(
            accu_mem_complex.begin(),
            miller::index<>(rfft.n_real()),
            miller::index<>(rfft.n_complex()));

        summation_eq15(space_group, miller_indices_f_obs,
          interm.m.const_ref(), f_part, fc_map, accu);
        rfft.backward(accu_mem_complex);
        maptbx::copy(accu_mem_real_const_ref, target_map_.ref());

        accu_mem_complex.fill(0);
        summation_eq14(space_group, miller_indices_f_obs,
          interm.m.const_ref(), f_part, fc_map, accu);
        rfft.backward(accu_mem_complex);
        combination_eq13(interm, accu_mem_real_const_ref, target_map_.ref());

        accu_mem_complex.fill(0);
        summation_eq14(space_group, miller_indices_f_obs,
          interm.m_d_i_obs.const_ref(), f_part, fc_map, accu);
        rfft.backward(accu_mem_complex);
        combination_eq12(interm, accu_mem_real_const_ref, target_map_.ref());
      }

      af::versa<FloatType, af::c_grid<3> >
      target_map() const { return target_map_; }

    protected:
      af::versa<FloatType, af::c_grid<3> > target_map_;
  };

}} // namespace cctbx::translation_search

#endif // CCTBX_TRANSLATION_SEARCH_FAST_NV1995_H
