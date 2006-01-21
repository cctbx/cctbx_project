#ifndef CCTBX_TRANSLATION_SEARCH_FAST_NV1995_H
#define CCTBX_TRANSLATION_SEARCH_FAST_NV1995_H

#include <cctbx/translation_search/fast_nv1995/combinations.h>
#include <cctbx/translation_search/fast_terms.h>

namespace cctbx { namespace translation_search {

  /*! \brief Fast computation of the translation function
      (correlation of intensities).
   */
  /*! Navaza, J. & Vernoslova, E. (1995). Acta Cryst. A51, 445-449.
   */
  template <typename FloatType=double>
  class fast_nv1995
  {
    public:
      //! Carries out the computation.
      fast_nv1995(
        af::int3 const& gridding,
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<miller::index<> > const& miller_indices_f_obs,
        af::const_ref<FloatType> const& f_obs,
        af::const_ref<std::complex<FloatType> > const& f_part,
        af::const_ref<miller::index<> > const& miller_indices_p1_f_calc,
        af::const_ref<std::complex<FloatType> > const& p1_f_calc)
      {
        fast_nv1995_detail::intermediates<FloatType> interm(
          space_group, anomalous_flag, miller_indices_f_obs, f_obs);
        fast_terms<FloatType> terms(
          gridding, anomalous_flag, miller_indices_p1_f_calc, p1_f_calc);
        target_map_ = terms.summation(
          space_group, miller_indices_f_obs,
          interm.m.const_ref(), f_part,
          true).fft().accu_real_copy();
        terms.summation(
          space_group, miller_indices_f_obs,
          interm.m.const_ref(), f_part,
          false).fft();
        fast_nv1995_detail::combination_eq13(
          interm, terms.accu_real_const_ref(), target_map_.ref());
        terms.summation(
          space_group, miller_indices_f_obs,
          interm.m_d_i_obs.const_ref(), f_part,
          false).fft();
        fast_nv1995_detail::combination_eq12(
          terms.accu_real_const_ref(), target_map_.ref());
      }

      //! Final result (correlation of intensities).
      af::versa<FloatType, af::c_grid<3> >
      target_map() const { return target_map_; }

    protected:
      af::versa<FloatType, af::c_grid<3> > target_map_;
  };

}} // namespace cctbx::translation_search

#endif // CCTBX_TRANSLATION_SEARCH_FAST_NV1995_H
