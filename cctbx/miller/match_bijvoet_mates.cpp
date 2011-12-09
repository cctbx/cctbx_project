#include <cctbx/miller/match_bijvoet_mates.h>
#include <cctbx/error.h>
#include <map>

namespace cctbx { namespace miller {

  void
  match_bijvoet_mates::match_(sgtbx::reciprocal_space::asu const& asu,
                              bool assert_is_unique_set_under_symmetry)
  {
    typedef std::map<index<>, std::size_t, fast_less_than<> > lookup_map_type;
    lookup_map_type lookup_map;
    for(std::size_t i=0;i<miller_indices_.size();i++) {
      lookup_map_type::const_iterator l = lookup_map.find(miller_indices_[i]);
      if (assert_is_unique_set_under_symmetry && l != lookup_map.end()) {
        throw CCTBX_ERROR("miller array is not a unique set under symmetry");
      }
      lookup_map[miller_indices_[i]] = i;
    }
    std::vector<bool> paired_already(miller_indices_.size(), false);
    for(std::size_t i=0;i<miller_indices_.size();i++) {
      if (paired_already[i]) continue;
      if (miller_indices_[i].is_zero()) {
        singles_[0].push_back(i);
      }
      else {
        int asu_which = asu.which(miller_indices_[i]);
        CCTBX_ASSERT(asu_which != 0);
        lookup_map_type::const_iterator
          l = lookup_map.find(-miller_indices_[i]);
        if (l == lookup_map.end()) {
          if (asu_which > 0) {
            singles_[0].push_back(i);
          }
          else {
            singles_[1].push_back(i);
          }
        }
        else {
          if (asu_which > 0) {
            pairs_.push_back(af::tiny<std::size_t, 2>(i, l->second));
          }
          else {
            pairs_.push_back(af::tiny<std::size_t, 2>(l->second, i));
          }
          paired_already[l->second] = true;
        }
      }
    }
  }

  af::shared<std::size_t>
  match_bijvoet_mates::singles(char plus_or_minus) const
  {
    std::size_t j = plus_or_minus_index_(plus_or_minus);
    return singles_[j];
  }

  void
  match_bijvoet_mates::size_assert_intrinsic() const
  {
    CCTBX_ASSERT(miller_indices_.size() == size_processed());
  }

  void
  match_bijvoet_mates::size_assert(std::size_t sz) const
  {
    size_assert_intrinsic();
    CCTBX_ASSERT(sz == size_processed());
  }

  af::shared<std::size_t>
  match_bijvoet_mates::pairs_hemisphere_selection(char plus_or_minus) const
  {
    std::size_t j = plus_or_minus_index_(plus_or_minus);
    af::const_ref<pair_type> pairs_ref = pairs_.const_ref();
    af::shared<std::size_t> result((af::reserve(pairs_ref.size())));
    for(std::size_t i=0;i<pairs_ref.size();i++) {
      result.push_back(pairs_[i][j]);
    }
    return result;
  }

  af::shared<index<> >
  match_bijvoet_mates::miller_indices_in_hemisphere(char plus_or_minus) const
  {
    std::size_t j = plus_or_minus_index_(plus_or_minus);
    af::shared<index<> > result((af::reserve(pairs_.size())));
    for(std::size_t i=0;i<pairs_.size();i++) {
      result.push_back(miller_indices_[pairs_[i][j]]);
    }
    return result;
  }

  std::size_t
  match_bijvoet_mates::plus_or_minus_index_(char plus_or_minus) const
  {
    CCTBX_ASSERT(plus_or_minus == '+' || plus_or_minus == '-');
    size_assert_intrinsic();
    if (plus_or_minus == '-') return 1;
    return 0;
  }

}} // namespace cctbx::miller
