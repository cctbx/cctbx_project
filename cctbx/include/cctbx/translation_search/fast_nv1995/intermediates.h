/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     Oct 2002: Fragments from phenix/fast_translation/front_end.h (rwgk)
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_TRANSLATION_SEARCH_FAST_NV1995_INTERMEDIATES_H
#define CCTBX_TRANSLATION_SEARCH_FAST_NV1995_INTERMEDIATES_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/shared_reductions.h>

namespace cctbx { namespace translation_search { namespace fast_nv1995_detail {

  template <typename FloatTypeSummation>
  struct intermediates
  {
    typedef af::shared<FloatTypeSummation> array_type;

    template <typename FloatTypeFobs,
              typename FloatTypeFpart>
    intermediates(sgtbx::space_group const& space_group,
                  bool anomalous_flag,
                  af::const_ref<miller::index<> > const& miller_indices,
                  af::const_ref<FloatTypeFobs> const& f_obs,
                  af::const_ref<std::complex<FloatTypeFpart> > const& f_part)
    {
      CCTBX_ASSERT(f_obs.size() == miller_indices.size());
      CCTBX_ASSERT(   f_part.size() == 0
                   || f_part.size() == miller_indices.size());
      array_type i_obs((af::reserve(miller_indices.size())));
      for(std::size_t i=0;i<miller_indices.size();i++) {
        i_obs.push_back(f_obs[i] * f_obs[i]);
      }
      array_type mh((af::reserve(miller_indices.size())));
      for(std::size_t i=0;i<miller_indices.size();i++) {
        mh.push_back(
          space_group.multiplicity(miller_indices[i], anomalous_flag));
      }
      sum_mh = af::sum(mh);
      d_i_obs = i_obs - af::sum(mh * i_obs) / sum_mh;
      sum_mh_d_i_obs_sq = af::sum(mh * d_i_obs * d_i_obs);
      sum_mh_f_sq = 0;
      sum_mh_f_sq_f_sq = 0;
      sum_mh_f_sq_d_i_obs = 0;
      if (f_part.size()) {
        array_type f_sq((af::reserve(f_part.size())));
        for(std::size_t i=0;i<f_part.size();i++) {
          f_sq.push_back(std::norm(f_part[i]));
        }
        sum_mh_f_sq = af::sum(mh * f_sq);
        sum_mh_f_sq_f_sq = af::sum(mh * f_sq * f_sq);
        sum_mh_f_sq_d_i_obs = af::sum(mh * f_sq * d_i_obs);
      }
    }

    FloatTypeSummation sum_mh;
    array_type         d_i_obs;
    FloatTypeSummation sum_mh_d_i_obs_sq;
    FloatTypeSummation sum_mh_f_sq;
    FloatTypeSummation sum_mh_f_sq_f_sq;
    FloatTypeSummation sum_mh_f_sq_d_i_obs;
  };

}}} // namespace cctbx::translation_search::fast_nv1995_detail

#endif // CCTBX_TRANSLATION_SEARCH_FAST_NV1995_INTERMEDIATES_H
