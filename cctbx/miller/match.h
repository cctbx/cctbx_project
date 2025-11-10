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
    struct additive_sigma {
      NumType operator()(NumType const& x, NumType const& y)
      {
        return std::sqrt(x*x + y*y);
      }
    };

    template <typename NumType>
    struct average {
      NumType operator()(NumType const& x, NumType const& y)
      {
        return (x + y) / NumType(2);
      }
    };

    template <typename Op>
    struct pair_op {

      pair_op(af::const_ref<pair_type> const& pairs)
      : pairs_(pairs)
      {}

      template <typename num_t>
      af::shared<num_t>
      operator()(
        af::const_ref<num_t> const& a0,
        af::const_ref<num_t> const& a1) const
      {
        af::shared<num_t> result((af::reserve(pairs_.size())));
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
