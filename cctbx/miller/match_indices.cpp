#include <cctbx/miller/match_indices.h>
#include <cctbx/error.h>
#include <map>
#include <vector>

namespace cctbx { namespace miller {


  match_indices::match_indices(
    af::shared<index<> > const& miller_indices_0)
  {
    singles_are_valid_ = false;
    pairs_are_valid_ = false;

    miller_indices_[0] = miller_indices_0;

    for(std::size_t i=0;i<miller_indices_[0].size();i++) {
      lookup_map_[miller_indices_[0][i]] = i;
    }
  }

  match_indices::match_indices(
    af::shared<index<> > const& miller_indices_0,
    af::shared<index<> > const& miller_indices_1)
  :
    miller_indices_(miller_indices_0, miller_indices_1),
    singles_are_valid_(true),
    pairs_are_valid_(true)
  {
    if (miller_indices_[0].id() == miller_indices_[1].id()) {
      // short-cut if same array
      pairs_.reserve(miller_indices_[0].size());
      for(std::size_t i=0;i<miller_indices_[0].size();i++) {
        pairs_.push_back(af::tiny<std::size_t, 2>(i, i));
      }
      return;
    }
    for(std::size_t i=0;i<miller_indices_[1].size();i++) {
      lookup_map_[miller_indices_[1][i]] = i;
    }
    std::vector<bool> miller_indices_1_flags(miller_indices_[1].size(), false);
    for(std::size_t i=0;i<miller_indices_[0].size();i++) {
      lookup_map_type::const_iterator
        l = lookup_map_.find(miller_indices_[0][i]);
      if (l == lookup_map_.end()) {
        singles_[0].push_back(i);
      }
      else {
        pairs_.push_back(af::tiny<std::size_t, 2>(i, l->second));
        miller_indices_1_flags[l->second] = true;
      }
    }
    for(std::size_t i=0;i<miller_indices_[1].size();i++) {
      if (!miller_indices_1_flags[i]) singles_[1].push_back(i);
    }
  }

  void
  match_indices::match_cached(
      af::shared<index<> > const& miller_indices_1)
  {
    /*
     * Match against the previously supplied miller_indices_0. This is faster
     * when calling repeatedly with different miller_indices_1 because the
     * lookup map is only constructed once.
     * */
    singles_are_valid_ = true;
    pairs_are_valid_ = true;

    miller_indices_[1] = miller_indices_1;
    pairs_.clear();
    singles_[0].clear();
    singles_[1].clear();
    af::shared<long> temp_pairs(lookup_map_.size(), -1);

    if (miller_indices_[0].id() == miller_indices_[1].id()) {
      // short-cut if same array
      pairs_.reserve(miller_indices_[0].size());
      for(std::size_t i=0;i<miller_indices_[0].size();i++) {
        pairs_.push_back(af::tiny<std::size_t, 2>(i, i));
      }
      return;
    }

    singles_[0].reserve(miller_indices_[0].size());
    singles_[1].reserve(miller_indices_[1].size());
    if (miller_indices_[0].size() < miller_indices_[1].size())
      pairs_.reserve(miller_indices_[0].size());
    else
      pairs_.reserve(miller_indices_[1].size());

    for(std::size_t i=0; i<miller_indices_[1].size(); i++) {
      lookup_map_type::const_iterator
        l = lookup_map_.find(miller_indices_[1][i]);
      if (l != lookup_map_.end())
        temp_pairs[l->second] = i;
      else
        singles_[1].push_back(i);
    }

    for (std::size_t i=0; i<temp_pairs.size(); ++i) {
      if (temp_pairs[i] != -1)
        pairs_.push_back(af::tiny<std::size_t, 2>(i, temp_pairs[i]));
      else
        singles_[0].push_back(i);
    }
  }

  void match_indices::match_cached_fast(
      af::shared<index<> > const& miller_indices_1)
  {
    /*
     * As match_indices::match_cached, but with some extra optimizations.
     * Only pairs are found, singles are ignored. The order of the results
     * is different (ordered as found in miller_indices_1, not _0). No member
     * functions except pairs() can be called if the matching was done this
     * way.
     * */
    singles_are_valid_ = false;
    pairs_are_valid_ = true;

    pairs_.clear();
    if (miller_indices_[0].size() < miller_indices_[1].size())
      pairs_.reserve(miller_indices_[0].size());
    else
      pairs_.reserve(miller_indices_[1].size());

    if (miller_indices_[0].id() == miller_indices_1.id()) {
      // short-cut if same array
      pairs_.reserve(miller_indices_[0].size());
      for(std::size_t i=0;i<miller_indices_[0].size();i++) {
        pairs_.push_back(af::tiny<std::size_t, 2>(i, i));
      }
      return;
    }

    for(std::size_t i=0;i<miller_indices_1.size();i++) {
      lookup_map_type::const_iterator
        l = lookup_map_.find(miller_indices_1[i]);
      if (l != lookup_map_.end())
        pairs_.push_back(af::tiny<std::size_t, 2>(l->second, i));
    }
  }

  void
  match_indices::size_assert_intrinsic() const
  {
    CCTBX_ASSERT(singles_are_valid_);
    CCTBX_ASSERT(pairs_are_valid_);
    CCTBX_ASSERT(miller_indices_[0].size() == size_processed(0));
    CCTBX_ASSERT(miller_indices_[1].size() == size_processed(1));
  }

  void
  match_indices::size_assert_1(std::size_t sz, std::size_t i) const
  {
    size_assert_intrinsic();
    CCTBX_ASSERT(sz == size_processed(i));
  }

  void
  match_indices::size_assert_2(std::size_t sz_0, std::size_t sz_1) const
  {
    size_assert_intrinsic();
    CCTBX_ASSERT(sz_0 == size_processed(0));
    CCTBX_ASSERT(sz_1 == size_processed(1));
  }

  af::shared<bool>
  match_indices::pair_selection(std::size_t i_array) const
  {
    CCTBX_ASSERT(i_array <= 1);
    size_assert_intrinsic();
    af::shared<bool> result(miller_indices_[i_array].size(), false);
    for(std::size_t i=0;i<pairs_.size();i++) {
      result[pairs_[i][i_array]] = true;
    }
    return result;
  }

  af::shared<bool>
  match_indices::single_selection(std::size_t i_array) const
  {
    CCTBX_ASSERT(i_array <= 1);
    size_assert_intrinsic();
    af::shared<bool> result(miller_indices_[i_array].size(), false);
    for(std::size_t i=0;i<singles_[i_array].size();i++) {
      result[singles_[i_array][i]] = true;
    }
    return result;
  }

  af::shared<index<> >
  match_indices::paired_miller_indices(std::size_t i_array) const
  {
    CCTBX_ASSERT(i_array <= 1);
    size_assert_intrinsic();
    af::shared<index<> > result((af::reserve(pairs_.size())));
    for(std::size_t i=0;i<pairs_.size();i++) {
      result.push_back(miller_indices_[i_array][pairs_[i][i_array]]);
    }
    return result;
  }

  af::shared<std::size_t>
  match_indices::permutation() const
  {
    size_assert_intrinsic();
    CCTBX_ASSERT(!singles_[0].size());
    af::shared<std::size_t> result((af::reserve(pairs_.size())));
    for(std::size_t i=0;i<pairs_.size();i++) {
      CCTBX_ASSERT(pairs_[i][0] == i);
      result.push_back(pairs_[i][1]);
    }
    return result;
  }

}} // namespace cctbx::miller
