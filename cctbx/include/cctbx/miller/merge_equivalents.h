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
      self.data.push_back(self.merge(
        current_index, &unmerged_data[group_begin], n));
      self.redundancies.push_back(n);
    }
  };

  namespace merge_equivalents {

    template <
      typename DerivedType,
      typename DataElementType,
      typename FloatType>
    void
    compute_r_fractors(
      DerivedType& self,
      const DataElementType* data_group,
      std::size_t n,
      FloatType const& result)
    {
      FloatType sum_num = scitbx::fn::absolute(data_group[0] - result);
      FloatType sum_den = scitbx::fn::absolute(data_group[0]);
      for(std::size_t i=1;i<n;i++) {
        sum_num += scitbx::fn::absolute(data_group[i] - result);
        sum_den += scitbx::fn::absolute(data_group[i]);
      }
      if (sum_den == 0) self.r_linear.push_back(0);
      else self.r_linear.push_back(sum_num / sum_den);
      //
      sum_num = scitbx::fn::pow2(data_group[0] - result);
      sum_den = scitbx::fn::pow2(data_group[0]);
      for(std::size_t i=1;i<n;i++) {
        sum_num += scitbx::fn::pow2(data_group[i] - result);
        sum_den += scitbx::fn::pow2(data_group[i]);
      }
      if (sum_den == 0) self.r_square.push_back(0);
      else self.r_square.push_back(sum_num / sum_den);
    }

  } // namespace merge_equivalents

  template <typename DataElementType, typename FloatType>
  struct merge_equivalents_generic : merge_equivalents_impl<DataElementType>
  {
    merge_equivalents_generic() {}

    merge_equivalents_generic(
      af::const_ref<index<> > const& unmerged_indices,
      af::const_ref<DataElementType> const& unmerged_data)
    {
      merge_equivalents_impl<DataElementType>
        ::loop_over_groups(*this, unmerged_indices, unmerged_data);
    }

    af::shared<index<> > indices;
    af::shared<DataElementType> data;
    af::shared<int> redundancies;

    DataElementType
    merge(
      miller::index<> const& /*current_index*/,
      const DataElementType* data_group, std::size_t n)
    {
      DataElementType result = data_group[0];
      for(std::size_t i=1;i<n;i++) result += data_group[i];
      return result / static_cast<FloatType>(n);
    }
  };

  template <typename IntegralType>
  struct merge_equivalents_exact : merge_equivalents_impl<IntegralType>
  {
    merge_equivalents_exact() {}

    merge_equivalents_exact(
      af::const_ref<index<> > const& unmerged_indices,
      af::const_ref<IntegralType> const& unmerged_data)
    {
      merge_equivalents_impl<IntegralType>
        ::loop_over_groups(*this, unmerged_indices, unmerged_data);
    }

    af::shared<index<> > indices;
    af::shared<IntegralType> data;
    af::shared<int> redundancies;

    IntegralType
    merge(
      miller::index<> const& current_index,
      const IntegralType* data_group, std::size_t n)
    {
      for(std::size_t i=1;i<n;i++) {
        if (data_group[i] != data_group[0]) {
          char buf[128];
          std::sprintf(buf,
            "merge_equivalents_exact:"
            " incompatible flags for hkl = (%d, %d, %d)",
            current_index[0], current_index[1], current_index[2]);
          throw error(buf);
        }
      }
      return data_group[0];
    }
  };

  template <typename FloatType=double>
  struct merge_equivalents_real : merge_equivalents_impl<FloatType>
  {
    merge_equivalents_real() {}

    merge_equivalents_real(
      af::const_ref<index<> > const& unmerged_indices,
      af::const_ref<FloatType> const& unmerged_data)
    {
      merge_equivalents_impl<FloatType>
        ::loop_over_groups(*this, unmerged_indices, unmerged_data);
    }

    af::shared<index<> > indices;
    af::shared<FloatType> data;
    af::shared<int> redundancies;
    //! r_linear = sum(abs(data - mean(data))) / sum(abs(data))
    af::shared<FloatType> r_linear;
    //! r_square = sum((data - mean(data))**2) / sum(data**2)
    af::shared<FloatType> r_square;

    FloatType
    merge(
      miller::index<> const& /*current_index*/,
      const FloatType* data_group, std::size_t n)
    {
      FloatType result = data_group[0];
      for(std::size_t i=1;i<n;i++) result += data_group[i];
      result /= static_cast<FloatType>(n);
      merge_equivalents::compute_r_fractors(*this, data_group, n, result);
      return result;
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
        af::const_ref<FloatType> const& unmerged_sigmas)
      {
        CCTBX_ASSERT(unmerged_data.size() == unmerged_indices.size());
        CCTBX_ASSERT(unmerged_sigmas.size() == unmerged_indices.size());
        init(unmerged_indices, unmerged_data, unmerged_sigmas);
      }

      af::shared<index<> > indices;
      af::shared<FloatType> data;
      af::shared<FloatType> sigmas;
      af::shared<int> redundancies;
      //! r_linear = sum(abs(data - mean(data))) / sum(abs(data))
      af::shared<FloatType> r_linear;
      //! r_square = sum((data - mean(data))**2) / sum(data**2)
      af::shared<FloatType> r_square;

    protected:
      void
      init(
        af::const_ref<index<> > const& unmerged_indices,
        af::const_ref<FloatType> const& unmerged_data,
        af::const_ref<FloatType> const& unmerged_sigmas)
      {
        if (unmerged_indices.size() == 0) return;
        std::vector<FloatType> values;
        std::vector<FloatType> weights;
        std::size_t group_begin = 0;
        std::size_t group_end = 1;
        for(;group_end<unmerged_indices.size();group_end++) {
          if (unmerged_indices[group_end] != unmerged_indices[group_begin]) {
            process_group(group_begin, group_end,
                          unmerged_indices[group_begin],
                          unmerged_data, unmerged_sigmas,
                          values, weights);
            group_begin = group_end;
          }
        }
        process_group(group_begin, group_end,
                      unmerged_indices[group_begin],
                      unmerged_data, unmerged_sigmas,
                      values, weights);
      }

      void
      process_group(std::size_t group_begin,
                    std::size_t group_end,
                    index<> const& current_index,
                    af::const_ref<FloatType> const& unmerged_data,
                    af::const_ref<FloatType> const& unmerged_sigmas,
                    std::vector<FloatType>& values,
                    std::vector<FloatType>& weights)
      {
        std::size_t n = group_end - group_begin;
        if (n == 0) return;
        indices.push_back(current_index);
        values.clear();
        values.reserve(n);
        weights.clear();
        weights.reserve(n);
        FloatType unmerged_sigma = 0;
        for(std::size_t i=0;i<n;i++) {
          FloatType ss = scitbx::fn::pow2(unmerged_sigmas[group_begin+i]);
          if (ss > 0) {
            values.push_back(unmerged_data[group_begin+i]);
            weights.push_back(1 / ss);
            unmerged_sigma = unmerged_sigmas[group_begin+i];
          }
        }
        if (values.size() == 0) {
          data.push_back(0);
          sigmas.push_back(0);
        }
        else if (values.size() == 1) {
          data.push_back(values[0]);
          sigmas.push_back(unmerged_sigma);
        }
        else {
          af::const_ref<FloatType> data_group(
            &*values.begin(), values.size());
          af::const_ref<FloatType> weights_group(
            &*weights.begin(), weights.size());
          scitbx::math::mean_and_variance<FloatType> mv(
            data_group, weights_group);
          data.push_back(mv.mean());
          sigmas.push_back(
            std::sqrt(
              std::max(
                mv.gsl_stats_wvariance()/values.size(),
                1/mv.sum_weights())));
        }
        redundancies.push_back(n);
        merge_equivalents::compute_r_fractors(
          *this, &unmerged_data[group_begin], n, data.back());
      }
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_MERGE_EQUIVALENTS_H
