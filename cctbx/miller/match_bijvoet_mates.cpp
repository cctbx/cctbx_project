/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Refactored parts of miller/miller_lib.cpp (rwgk)
 */

#include <cctbx/miller/match_bijvoet_mates.h>
#include <cctbx/error.h>
#include <map>

namespace cctbx { namespace miller {

  void
  match_bijvoet_mates::match_(sgtbx::reciprocal_space::asu const& asu)
  {
    typedef std::map<index<>, std::size_t> lookup_map_type;
    lookup_map_type lookup_map;
    for(std::size_t i=0;i<miller_indices_.size();i++) {
      lookup_map[miller_indices_[i]] = i;
    }
    std::vector<bool> paired_already(miller_indices_.size(), false);
    for(std::size_t i=0;i<miller_indices_.size();i++) {
      if (paired_already[i]) continue;
      lookup_map_type::const_iterator l = lookup_map.find(-miller_indices_[i]);
      if (l == lookup_map.end()) {
        singles_.push_back(i);
      }
      else {
        int asu_which = asu.which(miller_indices_[i]);
        CCTBX_ASSERT(asu_which != 0 || miller_indices_[i].is_zero());
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

  af::shared<bool>
  match_bijvoet_mates::hemisphere_selection(char plus_or_minus) const
  {
    CCTBX_ASSERT(plus_or_minus == '+' || plus_or_minus == '-');
    size_assert_intrinsic();
    std::size_t j = 0;
    if (plus_or_minus == '-') j = 1;
    af::shared<bool> result(miller_indices_.size(), false);
    for(std::size_t i=0;i<pairs_.size();i++) {
      result[pairs_[i][j]] = true;
    }
    return result;
  }

  af::shared<index<> >
  match_bijvoet_mates::miller_indices_in_hemisphere(char plus_or_minus) const
  {
    CCTBX_ASSERT(plus_or_minus == '+' || plus_or_minus == '-');
    size_assert_intrinsic();
    std::size_t j = 0;
    if (plus_or_minus == '-') j = 1;
    af::shared<index<> > result((af::reserve(pairs_.size())));
    for(std::size_t i=0;i<pairs_.size();i++) {
      result.push_back(miller_indices_[pairs_[i][j]]);
    }
    return result;
  }

}} // namespace cctbx::miller
