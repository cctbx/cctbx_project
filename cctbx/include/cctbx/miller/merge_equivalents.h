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

namespace cctbx { namespace miller {

  template <typename FloatType = double>
  class merge_equivalents
  {
    public:
      merge_equivalents() {}

      merge_equivalents(
        af::const_ref<index<> > const& unmerged_indices,
        af::const_ref<FloatType> const& unmerged_data,
        af::const_ref<FloatType> const& unmerged_sigmas)
      {
        CCTBX_ASSERT(unmerged_data.size() == unmerged_indices.size());
        CCTBX_ASSERT(unmerged_sigmas.size() == unmerged_indices.size());
        init(unmerged_indices, unmerged_data, unmerged_sigmas);
      }

      merge_equivalents(
        af::const_ref<index<> > const& unmerged_indices,
        af::const_ref<FloatType> const& unmerged_data)
      {
        CCTBX_ASSERT(unmerged_data.size() == unmerged_indices.size());
        af::const_ref<FloatType> unmerged_sigmas(0, 0);
        init(unmerged_indices, unmerged_data, unmerged_sigmas);
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
        af::const_ref<FloatType> const& unmerged_sigmas)
      {
        if (unmerged_indices.size() == 0) return;
        std::size_t group_begin = 0;
        std::size_t group_end = 1;
        for(;group_end<unmerged_indices.size();group_end++) {
          if (unmerged_indices[group_end] != unmerged_indices[group_begin]) {
            process_group(group_begin, group_end,
                          unmerged_indices[group_begin],
                          unmerged_data, unmerged_sigmas);
            group_begin = group_end;
          }
        }
        process_group(group_begin, group_end,
                      unmerged_indices[group_begin],
                      unmerged_data, unmerged_sigmas);
      }

      void
      process_group(std::size_t group_begin,
                    std::size_t group_end,
                    index<> const& current_index,
                    af::const_ref<FloatType> const& unmerged_data,
                    af::const_ref<FloatType> const& unmerged_sigmas)
      {
        std::size_t n = group_end - group_begin;
        if (n == 0) return;
        indices_.push_back(current_index);
        if (unmerged_sigmas.size() == 0) {
          FloatType sum_data = 0;
          for(std::size_t i=group_begin;i<group_end;i++) {
            sum_data += unmerged_data[i];
          }
          data_.push_back(sum_data / n);
        }
        else {
          FloatType sum_w_data = 0;
          FloatType sum_w = 0;
          for(std::size_t i=group_begin;i<group_end;i++) {
            CCTBX_ASSERT(unmerged_sigmas[i] != 0);
            FloatType s = unmerged_sigmas[i];
            FloatType w = 1 / (s*s);
            sum_w_data += unmerged_data[i] * w;
            sum_w += w;
          }
          CCTBX_ASSERT(sum_w > 0);
          data_.push_back(sum_w_data / sum_w);
          sigmas_.push_back(1 / std::sqrt(sum_w));
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
