#include <cctbx/miller/match_indices.h>
#include <cctbx/error.h>
#include <map>
#include <vector>

namespace cctbx { namespace miller {

  match_indices::match_indices(
    af::shared<index<> > const& miller_indices_0,
    af::shared<index<> > const& miller_indices_1)
  :
    miller_indices_(miller_indices_0, miller_indices_1)
  {
    if (miller_indices_[0].id() == miller_indices_[1].id()) {
      // short-cut if same array
      pairs_.reserve(miller_indices_[0].size());
      for(std::size_t i=0;i<miller_indices_[0].size();i++) {
        pairs_.push_back(af::tiny<std::size_t, 2>(i, i));
      }
      return;
    }
    typedef std::map<index<>, std::size_t, fast_less_than<> > lookup_map_type;
    lookup_map_type lookup_map;
    for(std::size_t i=0;i<miller_indices_[1].size();i++) {
      lookup_map[miller_indices_[1][i]] = i;
    }
    std::vector<bool> miller_indices_1_flags(miller_indices_[1].size(), false);
    for(std::size_t i=0;i<miller_indices_[0].size();i++) {
      lookup_map_type::const_iterator
        l = lookup_map.find(miller_indices_[0][i]);
      if (l == lookup_map.end()) {
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
  match_indices::size_assert_intrinsic() const
  {
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
    CCTBX_ASSERT(!have_singles());
    af::shared<std::size_t> result((af::reserve(pairs_.size())));
    for(std::size_t i=0;i<pairs_.size();i++) {
      CCTBX_ASSERT(pairs_[i][0] == i);
      result.push_back(pairs_[i][1]);
    }
    return result;
  }

}} // namespace cctbx::miller
