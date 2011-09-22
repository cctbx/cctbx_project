#ifndef CCTBX_MILLER_MATCH_BIJVOET_MATES_H
#define CCTBX_MILLER_MATCH_BIJVOET_MATES_H

#include <cctbx/miller/match.h>
#include <cctbx/sgtbx/reciprocal_space_asu.h>

namespace cctbx { namespace miller {

  class match_bijvoet_mates
  {
    public:
      match_bijvoet_mates() {}

      match_bijvoet_mates(
        sgtbx::space_group_type const& sg_type,
        af::shared<index<> > const& miller_indices,
        bool assert_is_unique_set_under_symmetry=true)
      :
        miller_indices_(miller_indices)
      {
        match_(sgtbx::reciprocal_space::asu(sg_type),
               assert_is_unique_set_under_symmetry);
      }

      match_bijvoet_mates(
        sgtbx::reciprocal_space::asu const& asu,
        af::shared<index<> > const& miller_indices,
        bool assert_is_unique_set_under_symmetry=true)
      :
        miller_indices_(miller_indices)
      {
        match_(asu, assert_is_unique_set_under_symmetry);
      }

      explicit
      match_bijvoet_mates(
        af::shared<index<> > const& miller_indices,
        bool assert_is_unique_set_under_symmetry=true)
      :
        miller_indices_(miller_indices)
      {
        match_(sgtbx::reciprocal_space::asu(sgtbx::space_group_type()),
               assert_is_unique_set_under_symmetry);
      }

      af::shared<pair_type>
      pairs() const
      {
        return pairs_;
      }

      af::shared<std::size_t>
      singles(char plus_or_minus) const;

      std::size_t
      n_singles() const
      {
        return singles_[0].size() + singles_[1].size();
      }

      /*! Not available in Python.
       */
      std::size_t
      size_processed() const
      {
        return 2 * pairs_.size() + n_singles();
      }

      /*! Not available in Python.
       */
      void
      size_assert_intrinsic() const;

      /*! Not available in Python.
       */
      void
      size_assert(std::size_t sz) const;

      af::shared<std::size_t>
      pairs_hemisphere_selection(char plus_or_minus) const;

      af::shared<std::size_t> const&
      singles_hemisphere_selection(char plus_or_minus) const
      {
        return singles_[plus_or_minus_index_(plus_or_minus)];
      }

      af::shared<index<> >
      miller_indices_in_hemisphere(char plus_or_minus) const;

      template <typename NumType>
      af::shared<NumType>
      minus(af::const_ref<NumType> const& data) const
      {
        size_assert(data.size());
        return detail::pair_op<std::minus<NumType> >
          (pairs_.const_ref())(data, data);
      }

      template <typename NumType>
      af::shared<NumType>
      additive_sigmas(af::const_ref<NumType> const& sigmas) const
      {
        size_assert(sigmas.size());
        return detail::pair_op<detail::additive_sigma<NumType> >
          (pairs_.const_ref())(sigmas, sigmas);
      }

      template <typename NumType>
      af::shared<NumType>
      average(af::const_ref<NumType> const& data) const
      {
        size_assert(data.size());
        return detail::pair_op<detail::average<NumType> >
          (pairs_.const_ref())(data, data);
      }

    protected:
      void
      match_(sgtbx::reciprocal_space::asu const& asu,
             bool assert_is_unique_set_under_symmetry=true);

      std::size_t
      plus_or_minus_index_(char plus_or_minus) const;

      af::shared<index<> > miller_indices_;
      af::shared<pair_type> pairs_;
      af::tiny<af::shared<std::size_t>, 2> singles_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_MATCH_BIJVOET_MATES_H
