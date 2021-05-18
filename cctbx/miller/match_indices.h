#ifndef CCTBX_MILLER_MATCH_INDICES_H
#define CCTBX_MILLER_MATCH_INDICES_H

#include <cctbx/miller/match.h>
#include <cctbx/miller.h>
#include <map>

namespace cctbx { namespace miller {

  typedef std::map<index<>, std::size_t, fast_less_than<> > lookup_map_type;

  class match_indices
  {
    public:

      match_indices() {}

      match_indices(af::shared<index<> > const& indices_0);

      match_indices(af::shared<index<> > const& indices_0,
                    af::shared<index<> > const& indices_1);

      void
      match_cached(af::shared<index<> > const& indices_1);

      void
      match_cached_fast(af::shared<index<> > const& indices_1);

      af::shared<pair_type>
      pairs() const
      {
        CCTBX_ASSERT(pairs_are_valid_);
        return pairs_;
      }

      af::shared<std::size_t>
      singles(std::size_t i) const
      {
        CCTBX_ASSERT(singles_are_valid_);
        if (i) return singles_[1];
        return singles_[0];
      }

      bool
      have_singles() const
      {
        CCTBX_ASSERT(singles_are_valid_);
        return singles_[0].size() || singles_[1].size();
      }

      /*! Not available in Python.
       */
      std::size_t
      size_processed(std::size_t i) const
      {
        CCTBX_ASSERT(singles_are_valid_);
        CCTBX_ASSERT(pairs_are_valid_);
        return pairs_.size() + singles_[i].size();
      }

      /*! Not available in Python.
       */
      void
      size_assert_intrinsic() const;

      /*! Not available in Python.
       */
      void
      size_assert_1(std::size_t sz, std::size_t i) const;

      /*! Not available in Python.
       */
      void
      size_assert_2(std::size_t sz_0, std::size_t sz_1) const;

      af::shared<bool>
      pair_selection(std::size_t i_array) const;

      af::shared<bool>
      single_selection(std::size_t i_array) const;

      af::shared<index<> >
      paired_miller_indices(std::size_t i_array) const;

      af::shared<std::size_t>
      permutation() const;

      template <typename NumType>
      af::shared<NumType>
      plus(
        af::const_ref<NumType> const& data_0,
        af::const_ref<NumType> const& data_1) const
      {
        size_assert_2(data_0.size(), data_1.size());
        return detail::pair_op<std::plus<NumType> >
          (pairs_.const_ref())(data_0, data_1);
      }

      template <typename NumType>
      af::shared<NumType>
      minus(
        af::const_ref<NumType> const& data_0,
        af::const_ref<NumType> const& data_1) const
      {
        size_assert_2(data_0.size(), data_1.size());
        return detail::pair_op<std::minus<NumType> >
          (pairs_.const_ref())(data_0, data_1);
      }

      template <typename NumType>
      af::shared<NumType>
      multiplies(
        af::const_ref<NumType> const& data_0,
        af::const_ref<NumType> const& data_1) const
      {
        size_assert_2(data_0.size(), data_1.size());
        return detail::pair_op<std::multiplies<NumType> >
          (pairs_.const_ref())(data_0, data_1);
      }

      template <typename NumType>
      af::shared<NumType>
      divides(
        af::const_ref<NumType> const& data_0,
        af::const_ref<NumType> const& data_1) const
      {
        size_assert_2(data_0.size(), data_1.size());
        return detail::pair_op<std::divides<NumType> >
          (pairs_.const_ref())(data_0, data_1);
      }

      template <typename NumType>
      af::shared<NumType>
      additive_sigmas(
        af::const_ref<NumType> const& sigmas_0,
        af::const_ref<NumType> const& sigmas_1) const
      {
        size_assert_2(sigmas_0.size(), sigmas_1.size());
        return detail::pair_op<detail::additive_sigma<NumType> >
          (pairs_.const_ref())(sigmas_0, sigmas_1);
      }

    protected:
      af::tiny<af::shared<index<> >, 2> miller_indices_;
      af::shared<pair_type> pairs_;
      af::tiny<af::shared<std::size_t>, 2> singles_;
      lookup_map_type lookup_map_;
      bool singles_are_valid_;
      bool pairs_are_valid_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_MATCH_INDICES_H
