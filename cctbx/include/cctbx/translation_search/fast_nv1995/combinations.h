#ifndef CCTBX_TRANSLATION_SEARCH_FAST_NV1995_COMBINATIONS_H
#define CCTBX_TRANSLATION_SEARCH_FAST_NV1995_COMBINATIONS_H

#include <cctbx/translation_search/fast_nv1995/intermediates.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>

namespace cctbx { namespace translation_search { namespace fast_nv1995_detail {

  // Navaza & Vernoslova (1995), p. 447, Eq. (13)
  template <typename FloatType>
  void
  combination_eq13(
    intermediates<FloatType> const& interm,
    af::const_ref<FloatType, af::c_grid_padded<3> > const&
      result_eq14,
    af::ref<FloatType, af::c_grid<3> > const&
      target_map) // on input result of equation 15
  {
    CCTBX_ASSERT(        result_eq14.accessor().focus()
                 .all_eq(target_map.accessor()));
    typedef FloatType f_t;
    typedef typename af::c_grid_padded<3>::index_type index_type;
    index_type n = result_eq14.accessor().focus();
    index_type i;
    std::size_t j = 0;
    for(i[0]=0;i[0]<n[0];i[0]++)
    for(i[1]=0;i[1]<n[1];i[1]++)
    for(i[2]=0;i[2]<n[2];i[2]++, j++) {
      f_t r14p = result_eq14(i);
      f_t r15p = target_map[j];
      f_t d = r15p - r14p * r14p / interm.sum_m;
      if (d > f_t(0)) {
        target_map[j] = std::sqrt(d) * std::sqrt(interm.sum_m_d_i_obs_sq);
      }
      else {
        target_map[j] = f_t(0);
      }
    }
  }

  // Navaza & Vernoslova (1995), p. 447, Eq. (12)
  template <typename FloatType>
  void
  combination_eq12(
    af::const_ref<FloatType, af::c_grid_padded<3> > const&
      result_eq14_with_i_obs,
    af::ref<FloatType, af::c_grid<3> > const&
      target_map, // on input result of equation 13
    FloatType const& big_correlation=1.e6)
  {
    CCTBX_ASSERT(        result_eq14_with_i_obs.accessor().focus()
                 .all_eq(target_map.accessor()));
    typedef FloatType f_t;
    typedef typename af::c_grid_padded<3>::index_type index_type;
    index_type n = result_eq14_with_i_obs.accessor().focus();
    index_type i;
    std::size_t j = 0;
    for(i[0]=0;i[0]<n[0];i[0]++)
    for(i[1]=0;i[1]<n[1];i[1]++)
    for(i[2]=0;i[2]<n[2];i[2]++, j++) {
      f_t r14ip = result_eq14_with_i_obs(i);
      if (scitbx::fn::absolute(r14ip / big_correlation) < target_map[j]) {
        target_map[j] = r14ip / target_map[j];
      }
      else {
        target_map[j] = f_t(0);
      }
    }
  }

}}} // namespace cctbx::translation_search::fast_nv1995_detail

#endif // CCTBX_TRANSLATION_SEARCH_FAST_NV1995_COMBINATIONS_H
