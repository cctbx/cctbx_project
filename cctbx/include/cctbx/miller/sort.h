/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created from fragments of cctbx/miller.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_SORT_H
#define CCTBX_MILLER_SORT_H

#include <cctbx/miller.h>
#include <cctbx/indexed_value.h>
#include <scitbx/array_family/shared.h>

namespace cctbx { namespace miller {

  namespace detail {

    template <typename DataType,
              typename SortCmpFunctor>
    void
    sort_in_place(
      af::shared<index<> > indices,
      af::shared<DataType> data,
      SortCmpFunctor sort_op)
    {
      CCTBX_ASSERT(indices.size() == data.size());
      typedef indexed_value<index<> , DataType, SortCmpFunctor> ivalue_type;
      af::shared<ivalue_type> ivalues((af::reserve(indices.size())));
      for(std::size_t i=0;i<indices.size();i++) {
        ivalues.push_back(ivalue_type(indices[i], data[i]));
      }
      std::sort(ivalues.begin(), ivalues.end());
      for(std::size_t i=0;i<indices.size();i++) {
        indices[i] = ivalues[i].index;
        data[i] = ivalues[i].value;
      }
    }

  } // namespace detail

  template <typename DataType>
  void
  sort_in_place(
    af::shared<index<> > indices,
    af::shared<DataType> data,
    bool reverse=false)
  {
    if (reverse) {
      detail::sort_in_place(indices, data, std::greater<DataType>());
    }
    else {
      detail::sort_in_place(indices, data, std::less<DataType>());
    }
  }

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_SORT_H
