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
#include <cctbx/sgtbx/miller_asu.h>
#include <cctbx/array_family/shared.h>

namespace cctbx { namespace miller {

  namespace detail {

    template <typename NumType>
    struct additive_sigma : std::binary_function<NumType, NumType, NumType>
    {
      typedef std::binary_function<NumType, NumType, NumType> base_type;

      typedef typename base_type::second_argument_type second_argument_type;
      typedef typename base_type::first_argument_type first_argument_type;
      typedef typename base_type::result_type result_type;

      NumType operator()(NumType const& x, NumType const& y)
      {
        return std::sqrt(x*x + y*y);
      }
    };

    template <typename NumType>
    struct average : std::binary_function<NumType, NumType, NumType>
    {
      typedef std::binary_function<NumType, NumType, NumType> base_type;

      typedef typename base_type::second_argument_type second_argument_type;
      typedef typename base_type::first_argument_type first_argument_type;
      typedef typename base_type::result_type result_type;

      NumType operator()(NumType const& x, NumType const& y)
      {
        return (x + y) / NumType(2);
      }
    };

    typedef af::tiny<std::size_t, 2> pair_type;

    template <typename Op>
    struct pair_op
    {
      typedef typename Op::second_argument_type second_argument_type;
      typedef typename Op::first_argument_type first_argument_type;
      typedef typename Op::result_type result_type;

      pair_op(af::shared<pair_type> pairs)
      : pairs_(pairs)
      {}

      af::shared<result_type>
      operator()(
        af::shared<first_argument_type> a0,
        af::shared<second_argument_type> a1) const
      {
        af::shared<result_type> result;
        result.reserve(pairs_.size());
        for(std::size_t i=0;i<pairs_.size();i++) {
          result.push_back(Op()(a0[pairs_[i][0]], a1[pairs_[i][1]]));
        }
        return result;
      }

      af::shared<pair_type> pairs_;
    };

  } // namespace detail

  class join_sets
  {
    public:
      join_sets() {}

      join_sets(af::shared<Index> miller_indices_0,
                af::shared<Index> miller_indices_1);

      af::shared<detail::pair_type> pairs() const
      {
        return pairs_;
      }

      af::shared<std::size_t> singles(std::size_t i) const
      {
        if (i) return singles_[1];
        return singles_[0];
      }

      bool have_singles() const
      {
        return singles_[0].size() || singles_[1].size();
      }

      std::size_t size_processed(std::size_t i) const
      {
        return pairs_.size() + singles_[i].size();
      }

      void size_assert_intrinsic() const
      {
        cctbx_assert(miller_indices_[0].size() == size_processed(0));
        cctbx_assert(miller_indices_[1].size() == size_processed(1));
      }

      void size_assert_1(std::size_t sz, std::size_t i) const
      {
        size_assert_intrinsic();
        cctbx_assert(sz == size_processed(i));
      }

      void size_assert_2(std::size_t sz_0, std::size_t sz_1) const
      {
        size_assert_intrinsic();
        cctbx_assert(sz_0 == size_processed(0));
        cctbx_assert(sz_1 == size_processed(1));
      }

      af::shared<Index>
      common_miller_indices() const;

      template <typename NumType>
      af::shared<NumType>
      plus(
        af::shared<NumType> data_0,
        af::shared<NumType> data_1) const
      {
        size_assert_2(data_0.size(), data_1.size());
        return detail::pair_op<std::plus<NumType> >(pairs_)(data_0, data_1);
      }

      template <typename NumType>
      af::shared<NumType>
      minus(
        af::shared<NumType> data_0,
        af::shared<NumType> data_1) const
      {
        size_assert_2(data_0.size(), data_1.size());
        return detail::pair_op<std::minus<NumType> >(pairs_)(data_0, data_1);
      }

      template <typename NumType>
      af::shared<NumType>
      multiplies(
        af::shared<NumType> data_0,
        af::shared<NumType> data_1) const
      {
        size_assert_2(data_0.size(), data_1.size());
        return
          detail::pair_op<std::multiplies<NumType> >(pairs_)(data_0, data_1);
      }

      template <typename NumType>
      af::shared<NumType>
      divides(
        af::shared<NumType> data_0,
        af::shared<NumType> data_1) const
      {
        size_assert_2(data_0.size(), data_1.size());
        return detail::pair_op<std::divides<NumType> >(pairs_)(data_0, data_1);
      }

      template <typename NumType>
      af::shared<NumType>
      additive_sigmas(
        af::shared<NumType> sigmas_0,
        af::shared<NumType> sigmas_1) const
      {
        size_assert_2(sigmas_0.size(), sigmas_1.size());
        return detail::pair_op<detail::additive_sigma<NumType> >(pairs_)(
          sigmas_0, sigmas_1);
      }

    protected:
      af::tiny<af::shared<Index>, 2> miller_indices_;
      af::shared<detail::pair_type> pairs_;
      af::tiny<af::shared<std::size_t>, 2> singles_;
  };

  class join_bijvoet_mates
  {
    public:
      join_bijvoet_mates() {}

      join_bijvoet_mates(
        sgtbx::SpaceGroupInfo const& sginfo,
        af::shared<Index> miller_indices)
      : miller_indices_(miller_indices)
      {
        join_(sgtbx::ReciprocalSpaceASU(sginfo));
      }

      join_bijvoet_mates(
        sgtbx::ReciprocalSpaceASU const& asu,
        af::shared<Index> miller_indices)
      : miller_indices_(miller_indices)
      {
        join_(asu);
      }

      explicit
      join_bijvoet_mates(
        af::shared<Index> miller_indices)
      : miller_indices_(miller_indices)
      {
        join_(sgtbx::ReciprocalSpaceASU(sgtbx::SpaceGroupInfo()));
      }

      af::shared<detail::pair_type> pairs() const
      {
        return pairs_;
      }

      af::shared<std::size_t> singles() const
      {
        return singles_;
      }

      bool have_singles() const {
        return singles_.size();
      }

      std::size_t size_processed() const
      {
        return 2 * pairs_.size() + singles_.size();
      }

      void size_assert_intrinsic() const
      {
        cctbx_assert(miller_indices_.size() == size_processed());
      }

      void size_assert(std::size_t sz) const
      {
        size_assert_intrinsic();
        cctbx_assert(sz == size_processed());
      }

      af::shared<Index>
      miller_indices_in_hemisphere(char plus_or_minus) const;

      template <typename NumType>
      af::shared<NumType>
      minus(af::shared<NumType> data) const
      {
        size_assert(data.size());
        return detail::pair_op<std::minus<NumType> >(pairs_)(data, data);
      }

      template <typename NumType>
      af::shared<NumType>
      additive_sigmas(af::shared<NumType> sigmas) const
      {
        size_assert(sigmas.size());
        return detail::pair_op<detail::additive_sigma<NumType> >(pairs_)(
          sigmas, sigmas);
      }

      template <typename NumType>
      af::shared<NumType>
      average(af::shared<NumType> data) const
      {
        size_assert(data.size());
        return detail::pair_op<detail::average<NumType> >(pairs_)(
          data, data);
      }

    protected:
      void join_(sgtbx::ReciprocalSpaceASU const& asu);

      af::shared<Index> miller_indices_;
      af::shared<detail::pair_type> pairs_;
      af::shared<std::size_t> singles_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_JOIN_H
