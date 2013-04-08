// -*- Mode: C++; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-

#include <cctbx/miller/match_multi_indices.h>
#include <cctbx/error.h>
#include <map>
#include <vector>

namespace cctbx { namespace miller {

  match_multi_indices::match_multi_indices(
    af::shared<index<> > const& miller_indices_set,
    af::shared<index<> > const& miller_indices_multiset)
  :
    miller_indices_(miller_indices_set, miller_indices_multiset),
    number_of_matches_(af::shared<std::size_t>(miller_indices_set.size()),
                       af::shared<std::size_t>(miller_indices_multiset.size()))
  {
    typedef std::map<index<>, std::size_t, fast_less_than<> > lookup_map_type;
    lookup_map_type lookup_map;
    for(std::size_t i=0;i<miller_indices_[0].size();i++) {
      std::pair<lookup_map_type::const_iterator, bool>
        s = lookup_map.insert(std::make_pair(miller_indices_[0][i], i));
      CCTBX_ASSERT(s.second);
    }
    for(std::size_t i=0;i<miller_indices_[1].size();i++) {
      lookup_map_type::const_iterator
        l = lookup_map.find(miller_indices_[1][i]);
      if(l!=lookup_map.end()) {
        for(;l->first==miller_indices_[1][i] && l!=lookup_map.end();l++) {
          pairs_.push_back(af::tiny<std::size_t, 2>(l->second, i));
          number_of_matches_[0][l->second]++;
          number_of_matches_[1][i]++;
        }
      }
    }
  }

  bool
  match_multi_indices::have_singles() const
  {
    for(std::size_t i=0;i<2;i++) {
      for(std::size_t j=0;j<number_of_matches_[i].size();j++) {
        if(number_of_matches_[i][j]==0)
          return true;
      }
    }
    return false;
  }

  af::shared<std::size_t>
  match_multi_indices::singles(std::size_t i_array) const
  {
    if(i_array!=0)
      i_array = 1;
    af::shared<std::size_t>
      result((af::reserve(number_of_matches_[i_array].size())));
    for(std::size_t j=0;j<number_of_matches_[i_array].size();j++) {
      if(number_of_matches_[i_array][j]==0)
        result.push_back(j);
    }
    return result;
  }

  af::shared<bool>
  match_multi_indices::pair_selection(std::size_t i_array) const
  {
    CCTBX_ASSERT(i_array <= 1);
    af::shared<bool> result(miller_indices_[i_array].size(), false);
    for(std::size_t i=0;i<miller_indices_[i_array].size();i++) {
      result[i] = number_of_matches_[i_array][i] > 0;
    }
    return result;
  }

  af::shared<bool>
  match_multi_indices::single_selection(std::size_t i_array) const
  {
    CCTBX_ASSERT(i_array <= 1);
    af::shared<bool> result(miller_indices_[i_array].size(), false);
    for(std::size_t i=0;i<miller_indices_[i_array].size();i++) {
      result[i] = number_of_matches_[i_array][i] == 0;
    }
    return result;
  }

  af::shared<index<> >
  match_multi_indices::paired_miller_indices(std::size_t i_array) const
  {
    CCTBX_ASSERT(i_array <= 1);
    af::shared<index<> > result((af::reserve(pairs_.size())));
    for(std::size_t i=0;i<pairs_.size();i++) {
      if(number_of_matches_[i_array][i] > 0) {
        result.push_back(miller_indices_[i_array][i]);
      }
    }
    return result;
  }

}} // namespace cctbx::miller
