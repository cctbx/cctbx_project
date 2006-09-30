#ifndef CCTBX_MILLER_MATCH_H
#define CCTBX_MILLER_MATCH_H

#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>
#include <functional>
#include <cmath>

namespace cctbx { namespace miller {

  typedef af::tiny<std::size_t, 2> pair_type;

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

    template <typename Op>
    struct pair_op
    {
      typedef typename Op::second_argument_type second_argument_type;
      typedef typename Op::first_argument_type first_argument_type;
      typedef typename Op::result_type result_type;

      pair_op(af::const_ref<pair_type> const& pairs)
      : pairs_(pairs)
      {}

      af::shared<result_type>
      operator()(
        af::const_ref<first_argument_type> const& a0,
        af::const_ref<second_argument_type> const& a1) const
      {
        af::shared<result_type> result((af::reserve(pairs_.size())));
        for(std::size_t i=0;i<pairs_.size();i++) {
          result.push_back(Op()(a0[pairs_[i][0]], a1[pairs_[i][1]]));
        }
        return result;
      }

      af::const_ref<pair_type> pairs_;
    };

  } // namespace detail

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_MATCH_H
