/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Jan: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_SORT_H
#define SCITBX_ARRAY_FAMILY_SORT_H

#include <scitbx/indexed_value.h>
#include <scitbx/array_family/shared.h>

namespace scitbx { namespace af {

  namespace detail {

    template <typename DataType,
              typename SortCmpFunctor>
    shared<std::size_t>
    sort_permutation(
      const_ref<DataType> const& data,
      SortCmpFunctor sort_op)
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

  } // namespace detail

  template <typename DataType>
  shared<std::size_t>
  sort_permutation(
    const_ref<DataType> const& data,
    bool reverse=false)
  {
    if (reverse) {
      return detail::sort_permutation(data, std::greater<DataType>());
    }
    else {
      return detail::sort_permutation(data, std::less<DataType>());
    }
  }

  template <typename DataType>
  shared<DataType>
  shuffle(
    const_ref<DataType> const& data,
    const_ref<std::size_t> const& permutation)
  {
    shared<DataType> result((reserve(data.size())));
    for(std::size_t i=0;i<permutation.size();i++) {
      std::size_t permutation_i = permutation[i];
      SCITBX_ASSERT(permutation_i < data.size());
      result.push_back(data[permutation_i]);
    }
    return result;
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_SORT_H
