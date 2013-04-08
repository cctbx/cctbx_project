// -*- Mode: C++; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-

#ifndef CCTBX_MILLER_MATCH_MULTI_INDICES_H
#define CCTBX_MILLER_MATCH_MULTI_INDICES_H

#include <cctbx/miller/match.h>
#include <cctbx/miller.h>

namespace cctbx { namespace miller {

  class match_multi_indices
  {
    public:
      match_multi_indices() {}

      match_multi_indices(af::shared<index<> > const& miller_indices_unique,
                          af::shared<index<> > const& miller_indices);

      af::shared<std::size_t>
      number_of_matches(std::size_t i_array) const
      {
        CCTBX_ASSERT(i_array <= 1);
        return number_of_matches_[i_array];
      }

      bool
      have_singles() const;

      af::shared<pair_type>
      pairs() const
      {
        return pairs_;
      }

      af::shared<std::size_t>
      singles(std::size_t i_array) const;

      af::shared<bool>
      pair_selection(std::size_t i_array) const;

      af::shared<bool>
      single_selection(std::size_t i_array) const;

      af::shared<index<> >
      paired_miller_indices(std::size_t i_array) const;

    protected:
      af::tiny<af::shared<index<> >, 2> miller_indices_;
      af::tiny<af::shared<std::size_t>, 2> number_of_matches_;
      af::shared<pair_type> pairs_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_MATCH_MULTI_INDICES_H
