/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created from fragments of cctbx/miller.h (R.W. Grosse-Kunstleve)
 */

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
        af::shared<index<> > const& miller_indices)
      :
        miller_indices_(miller_indices)
      {
        match_(sgtbx::reciprocal_space::asu(sg_type));
      }

      match_bijvoet_mates(
        sgtbx::reciprocal_space::asu const& asu,
        af::shared<index<> > const& miller_indices)
      :
        miller_indices_(miller_indices)
      {
        match_(asu);
      }

      explicit
      match_bijvoet_mates(
        af::shared<index<> > const& miller_indices)
      :
        miller_indices_(miller_indices)
      {
        match_(sgtbx::reciprocal_space::asu(sgtbx::space_group_type()));
      }

      af::shared<pair_type>
      pairs() const
      {
        return pairs_;
      }

      af::shared<std::size_t>
      singles() const
      {
        return singles_;
      }

      bool
      have_singles() const
      {
        return singles_.size();
      }

      /*! Not available in Python.
       */
      std::size_t
      size_processed() const
      {
        return 2 * pairs_.size() + singles_.size();
      }

      /*! Not available in Python.
       */
      void
      size_assert_intrinsic() const;

      /*! Not available in Python.
       */
      void
      size_assert(std::size_t sz) const;

      af::shared<bool>
      hemisphere_selection(char plus_or_minus) const;

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
      match_(sgtbx::reciprocal_space::asu const& asu);

      af::shared<index<> > miller_indices_;
      af::shared<pair_type> pairs_;
      af::shared<std::size_t> singles_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_MATCH_BIJVOET_MATES_H
