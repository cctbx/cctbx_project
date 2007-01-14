#ifndef SCITBX_ARRAY_FAMILY_RANGE_H
#define SCITBX_ARRAY_FAMILY_RANGE_H

#include <scitbx/array_family/shared.h>
#include <stdexcept>

namespace scitbx { namespace af {

  namespace range_args {

    template <typename IntType>
    struct no_check
    {
      static void check_start(IntType const&) {}
      static void check_stop(IntType const&) {}
    };

    template <typename IntType>
    struct unsigned_check
    {
      static void check_start(IntType const& start)
      {
        if (start < 0) throw std::runtime_error(
          "range start argument must not be negative.");
      }

      static void check_stop(IntType const& stop)
      {
        if (stop < 0) throw std::runtime_error(
          "range stop argument must not be negative.");
      }
    };

  } // namespace range_args

  template <typename ValueType,
            typename IntType,
            template<typename> class CheckType=range_args::no_check>
  struct range
  {
    // Python-2.5/Python/bltinmodule.c get_len_of_range()
    static
    std::size_t
    size(
      IntType const& low,
      IntType const& high,
      IntType const& step)
    {
      if (step == 0) {
        throw std::runtime_error("range step argument must not be zero.");
      }
      if (low < high) {
        return (  static_cast<std::size_t>(high)
                - static_cast<std::size_t>(low) - 1)
              / static_cast<std::size_t>(step) + 1;
      }
      return 0;
    }

    // Python-2.5/Python/bltinmodule.c builtin_range()
    static
    shared<ValueType>
    array(
      IntType const& start,
      IntType const& stop,
      IntType const& step)
    {
      CheckType<IntType>::check_start(start);
      CheckType<IntType>::check_stop(stop);
      shared<ValueType> result;
      std::size_t n;
      if (step < 0) n = size(stop, start, -step);
      else          n = size(start, stop, step);
      result.reserve(n);
      IntType ival = start;
      for (std::size_t i=0;i<n;i++) {
        result.push_back(static_cast<ValueType>(ival));
        ival += step;
      }
      return result;
    }

    static
    shared<ValueType>
    array(
      IntType const& start,
      IntType const& stop)
    {
      return array(start, stop, 1);
    }

    static
    shared<ValueType>
    array(
      IntType const& stop)
    {
      return array(0, stop, 1);
    }
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_RANGE_H
