// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created from fragments of cctbx/miller.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_JOIN_H
#define CCTBX_MILLER_JOIN_H

#include <vector>
#include <map>
#include <algorithm>
#include <cctbx/miller.h>
#include <cctbx/sgtbx/miller_asu.h>
#include <cctbx/array_family/shared.h>

namespace cctbx { namespace miller {

  class join_sets
  {
    public:
      join_sets() {}

      join_sets(af::shared<Index> a1,
                af::shared<Index> a2);

      af::shared<af::tiny<std::size_t, 2> > pairs() const {
        return pairs_;
      }

      af::shared<std::size_t> singles(std::size_t i) const {
        if (i) return singles_[1];
        return singles_[0];
      }

      bool have_singles() const {
        return singles_[0].size() || singles_[1].size();
      }

    protected:
      af::shared<af::tiny<std::size_t, 2> > pairs_;
      af::shared<std::size_t> singles_[2];
  };

  class join_bijvoet_mates
  {
    public:
      join_bijvoet_mates() {}

      join_bijvoet_mates(
        sgtbx::SpaceGroupInfo const& sginfo,
        af::shared<Index> miller_indices)
      {
        join_(sgtbx::ReciprocalSpaceASU(sginfo), miller_indices);
      }

      join_bijvoet_mates(
        sgtbx::ReciprocalSpaceASU const& asu,
        af::shared<Index> miller_indices)
      {
        join_(asu, miller_indices);
      }

      join_bijvoet_mates(
        af::shared<Index> miller_indices)
      {
        join_(
          sgtbx::ReciprocalSpaceASU(sgtbx::SpaceGroupInfo()), miller_indices);
      }

      af::shared<af::tiny<std::size_t, 2> > pairs() const {
        return pairs_;
      }

      af::shared<std::size_t> singles() const {
        return singles_;
      }

      bool have_singles() const {
        return singles_.size();
      }

      af::shared<Index>
      select(af::shared<Index> miller_indices, bool plus) const;

    protected:
      void join_(
        sgtbx::ReciprocalSpaceASU const& asu,
        af::shared<Index> miller_indices);

      af::shared<af::tiny<std::size_t, 2> > pairs_;
      af::shared<std::size_t> singles_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_JOIN_H
