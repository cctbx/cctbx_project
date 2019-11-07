#ifndef SCITBX_ARRAY_FAMILY_SORT_H
#define SCITBX_ARRAY_FAMILY_SORT_H

#include <scitbx/indexed_value.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/error.h>

namespace scitbx { namespace af {

  namespace detail {

    template <typename DataType,
              typename SortCmpFunctor>
    shared<std::size_t>
    sort_permutation(
      const_ref<DataType> const& data,
      SortCmpFunctor /*sort_op*/)
    {
      typedef indexed_value<
        std::size_t , DataType, SortCmpFunctor> ivalue_type;
      shared<std::size_t> result((reserve(data.size())));
      shared<ivalue_type> ivalues((reserve(data.size())));
      for(std::size_t i=0;i<data.size();i++) {
        ivalues.push_back(ivalue_type(i, data[i]));
      }
      std::sort(ivalues.begin(), ivalues.end());
      for(std::size_t i=0;i<data.size();i++) {
        result.push_back(ivalues[i].index);
      }
      return result;
    }

    template <typename DataType,
              typename SortCmpFunctor>
    shared<std::size_t>
    stable_sort_permutation(
      const_ref<DataType> const& data,
      SortCmpFunctor /*sort_op*/)
    {

      typedef indexed_value<
        std::size_t , DataType, SortCmpFunctor> ivalue_type;
      shared<std::size_t> result((reserve(data.size())));
      shared<ivalue_type> ivalues((reserve(data.size())));
      for(std::size_t i=0;i<data.size();i++) {
        ivalues.push_back(ivalue_type(i, data[i]));
      }
      std::stable_sort(ivalues.begin(), ivalues.end());
      for(std::size_t i=0;i<data.size();i++) {
        result.push_back(ivalues[i].index);
      }
      return result;
    }

  } // namespace detail

  template <typename DataType>
  shared<std::size_t>
  sort_permutation(
    const_ref<DataType> const& data,
    bool reverse=false,
    bool stable=true
    )
  {
    if (stable) {
      if (reverse) {
        return detail::stable_sort_permutation(data, std::greater<DataType>());
      }
      else {
        return detail::stable_sort_permutation(data, std::less<DataType>());
      }
    }
    else {
      if (reverse) {
        return detail::sort_permutation(data, std::greater<DataType>());
      }
      else {
        return detail::sort_permutation(data, std::less<DataType>());
      }
    }
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_SORT_H
