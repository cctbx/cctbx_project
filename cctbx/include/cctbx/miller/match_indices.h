/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created from fragments of cctbx/miller.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_MATCH_INDICES_H
#define CCTBX_MILLER_MATCH_INDICES_H

#include <cctbx/miller/match.h>
#include <cctbx/miller.h>

namespace cctbx { namespace miller {

  class match_indices
  {
    public:
      match_indices() {}

      match_indices(af::shared<index<> > const& indices_0,
                    af::shared<index<> > const& indices_1);

      af::shared<pair_type>
      pairs() const
      {
        return pairs_;
      }

      af::shared<std::size_t>
      singles(std::size_t i) const
      {
        if (i) return singles_[1];
        return singles_[0];
      }

      bool
      have_singles() const
      {
        return singles_[0].size() || singles_[1].size();
      }

      /*! Not available in Python.
       */
      std::size_t
      size_processed(std::size_t i) const
      {
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
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_MATCH_INDICES_H
