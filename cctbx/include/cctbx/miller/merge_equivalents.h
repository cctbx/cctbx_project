/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Apr: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_MERGE_EQUIVALENTS_H
#define CCTBX_MILLER_MERGE_EQUIVALENTS_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/math/mean_and_variance.h>

namespace cctbx { namespace miller {

  template <typename FloatType = double>
  class merge_equivalents
  {
    public:
      merge_equivalents() {}

      merge_equivalents(
        af::const_ref<index<> > const& unmerged_indices,
        af::const_ref<FloatType> const& unmerged_data,
        af::const_ref<FloatType> const& unmerged_weights)
      {
        CCTBX_ASSERT(unmerged_data.size() == unmerged_indices.size());
        CCTBX_ASSERT(unmerged_weights.size() == unmerged_indices.size());
        init(unmerged_indices, unmerged_data, unmerged_weights);
      }

      merge_equivalents(
        af::const_ref<index<> > const& unmerged_indices,
        af::const_ref<FloatType> const& unmerged_data)
      {
        CCTBX_ASSERT(unmerged_data.size() == unmerged_indices.size());
        af::const_ref<FloatType> unmerged_weights(0, 0);
        init(unmerged_indices, unmerged_data, unmerged_weights);
      }

      af::shared<index<> >
      indices() const { return indices_; }

      af::shared<FloatType>
      data() const { return data_; }

      af::shared<FloatType>
      sigmas() const { return sigmas_; }

      af::shared<int>
      redundancies() const { return redundancies_; }

    protected:
      void
      init(
        af::const_ref<index<> > const& unmerged_indices,
        af::const_ref<FloatType> const& unmerged_data,
        af::const_ref<FloatType> const& unmerged_weights)
      {
        if (unmerged_indices.size() == 0) return;
        std::size_t group_begin = 0;
        std::size_t group_end = 1;
        for(;group_end<unmerged_indices.size();group_end++) {
          if (unmerged_indices[group_end] != unmerged_indices[group_begin]) {
            process_group(group_begin, group_end,
                          unmerged_indices[group_begin],
                          unmerged_data, unmerged_weights);
            group_begin = group_end;
          }
        }
        process_group(group_begin, group_end,
                      unmerged_indices[group_begin],
                      unmerged_data, unmerged_weights);
      }

      void
      process_group(std::size_t group_begin,
                    std::size_t group_end,
                    index<> const& current_index,
                    af::const_ref<FloatType> const& unmerged_data,
                    af::const_ref<FloatType> const& unmerged_weights)
      {
        std::size_t n = group_end - group_begin;
        if (n == 0) return;
        indices_.push_back(current_index);
        if (n == 1) {
          data_.push_back(unmerged_data[group_begin]);
          if (unmerged_weights.size() != 0) {
            sigmas_.push_back(1/std::sqrt(unmerged_weights[group_begin]));
          }
        }
        else {
          af::const_ref<FloatType> data_group(
            &unmerged_data[group_begin], n);
          if (unmerged_weights.size() == 0) {
            scitbx::math::mean_and_variance<FloatType> mv(data_group);
            data_.push_back(mv.mean());
          }
          else {
            af::const_ref<FloatType> weights_group(
              &unmerged_weights[group_begin],n);
            scitbx::math::mean_and_variance<FloatType> mv(
              data_group, weights_group);
            data_.push_back(mv.mean());
            sigmas_.push_back(mv.standard_deviation());
          }
        }
        redundancies_.push_back(n);
      }

      af::shared<index<> > indices_;
      af::shared<FloatType> data_;
      af::shared<FloatType> sigmas_;
      af::shared<int> redundancies_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_MERGE_EQUIVALENTS_H
