#ifndef CCTBX_TRANSLATION_SEARCH_FAST_NV1995_INTERMEDIATES_H
#define CCTBX_TRANSLATION_SEARCH_FAST_NV1995_INTERMEDIATES_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/shared_reductions.h>

namespace cctbx { namespace translation_search { namespace fast_nv1995_detail {

  template <typename FloatType>
  struct intermediates
  {
    typedef af::shared<FloatType> array_type;

    intermediates(sgtbx::space_group const& space_group,
                  bool anomalous_flag,
                  af::const_ref<miller::index<> > const& miller_indices,
                  af::const_ref<FloatType> const& f_obs)
    {
      CCTBX_ASSERT(f_obs.size() == miller_indices.size());
      array_type i_obs((af::reserve(miller_indices.size())));
      for(std::size_t i=0;i<miller_indices.size();i++) {
        i_obs.push_back(f_obs[i] * f_obs[i]);
      }
      m.reserve(miller_indices.size());
      for(std::size_t i=0;i<miller_indices.size();i++) {
        m.push_back(
          space_group.multiplicity(miller_indices[i], anomalous_flag));
      }
      sum_m = af::sum(m);
      array_type d_i_obs = i_obs - af::sum(m * i_obs) / sum_m;
      m_d_i_obs = m * d_i_obs;
      sum_m_d_i_obs_sq = af::sum(m_d_i_obs * d_i_obs);
    }

    array_type m;
    FloatType sum_m;
    array_type m_d_i_obs;
    FloatType sum_m_d_i_obs_sq;
  };

}}} // namespace cctbx::translation_search::fast_nv1995_detail

#endif // CCTBX_TRANSLATION_SEARCH_FAST_NV1995_INTERMEDIATES_H
