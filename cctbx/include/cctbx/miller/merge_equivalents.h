#ifndef CCTBX_MILLER_MERGE_EQUIVALENTS_H
#define CCTBX_MILLER_MERGE_EQUIVALENTS_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/math/mean_and_variance.h>

namespace cctbx { namespace miller {

  template <typename DataElementType>
  struct merge_equivalents_impl
  {
    template <typename DerivedType>
    void
    loop_over_groups(
      DerivedType& self,
      af::const_ref<index<> > const& unmerged_indices,
      af::const_ref<DataElementType> const& unmerged_data)
    {
      CCTBX_ASSERT(unmerged_data.size() == unmerged_indices.size());
      if (unmerged_indices.size() == 0) return;
      std::size_t group_begin = 0;
      std::size_t group_end = 1;
      for(;group_end<unmerged_indices.size();group_end++) {
        if (unmerged_indices[group_end] != unmerged_indices[group_begin]) {
          process_group(
            self, group_begin, group_end,
            unmerged_indices[group_begin], unmerged_data);
          group_begin = group_end;
        }
      }
      process_group(
        self, group_begin, group_end,
        unmerged_indices[group_begin], unmerged_data);
    }

    template <typename DerivedType>
    void
    process_group(
      DerivedType& self,
      std::size_t group_begin,
      std::size_t group_end,
      index<> const& current_index,
      af::const_ref<DataElementType> const& unmerged_data)
    {
      std::size_t n = group_end - group_begin;
      if (n == 0) return;
      self.indices.push_back(current_index);
      if (n == 1) {
        self.data.push_back(unmerged_data[group_begin]);
      }
      else {
        self.data.push_back(self.merge(
          current_index, &unmerged_data[group_begin], n));
      }
      self.redundancies.push_back(n);
    }
  };

  template <typename DataType, typename FloatType>
  struct merge_equivalents_generic : merge_equivalents_impl<DataType>
  {
    merge_equivalents_generic() {}

    merge_equivalents_generic(
      af::const_ref<index<> > const& unmerged_indices,
      af::const_ref<DataType> const& unmerged_data)
    {
      merge_equivalents_impl<DataType>
        ::loop_over_groups(*this, unmerged_indices, unmerged_data);
    }

    af::shared<index<> > indices;
    af::shared<DataType> data;
    af::shared<int> redundancies;

    DataType
    merge(
      miller::index<> const& current_index,
      const DataType* data_group, std::size_t n)
    {
      DataType result = data_group[0];
      for(std::size_t i=1;i<n;i++) result += data_group[i];
      return result / static_cast<FloatType>(n);
    }
  };

  struct merge_equivalents_bool : merge_equivalents_impl<bool>
  {
    merge_equivalents_bool() {}

    merge_equivalents_bool(
      af::const_ref<index<> > const& unmerged_indices,
      af::const_ref<bool> const& unmerged_data)
    {
      merge_equivalents_impl<bool>
        ::loop_over_groups(*this, unmerged_indices, unmerged_data);
    }

    af::shared<index<> > indices;
    af::shared<bool> data;
    af::shared<int> redundancies;

    bool
    merge(
      miller::index<> const& current_index,
      const bool* data_group, std::size_t n)
    {
      for(std::size_t i=1;i<n;i++) {
        if (data_group[i] != data_group[0]) {
          char buf[128];
          std::sprintf(buf,
            "merge_equivalents_bool:"
            " incompatible flags for hkl = (%d, %d, %d)",
            current_index[0], current_index[1], current_index[2]);
          throw error(buf);
        }
      }
      return data_group[0];
    }
  };

  template <typename FloatType=double>
  class merge_equivalents_obs
  {
    public:
      merge_equivalents_obs() {}

      merge_equivalents_obs(
        af::const_ref<index<> > const& unmerged_indices,
        af::const_ref<FloatType> const& unmerged_data,
        af::const_ref<FloatType> const& unmerged_weights)
      {
        CCTBX_ASSERT(unmerged_data.size() == unmerged_indices.size());
        CCTBX_ASSERT(unmerged_weights.size() == unmerged_indices.size());
        init(unmerged_indices, unmerged_data, unmerged_weights);
      }

      af::shared<index<> > indices;
      af::shared<FloatType> data;
      af::shared<FloatType> sigmas;
      af::shared<int> redundancies;

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
        indices.push_back(current_index);
        if (n == 1) {
          data.push_back(unmerged_data[group_begin]);
          sigmas.push_back(1/std::sqrt(unmerged_weights[group_begin]));
        }
        else {
          af::const_ref<FloatType> data_group(
            &unmerged_data[group_begin], n);
          af::const_ref<FloatType> weights_group(
            &unmerged_weights[group_begin],n);
          scitbx::math::mean_and_variance<FloatType> mv(
            data_group, weights_group);
          data.push_back(mv.mean());
          sigmas.push_back(mv.conservative_standard_deviation());
        }
        redundancies.push_back(n);
      }
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_MERGE_EQUIVALENTS_H
