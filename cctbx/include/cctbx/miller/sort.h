// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created from fragments of cctbx/miller.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_SORT_H
#define CCTBX_MILLER_SORT_H

#include <algorithm>
#include <cctbx/miller.h>
#include <cctbx/indexed_value.h>
#include <cctbx/array_family/shared.h>

namespace cctbx { namespace miller {

  template <typename DataType,
            typename SortCmpFunctor>
  void
  inplace_sort(
    af::shared<Index> miller_indices,
    af::shared<DataType> data,
    SortCmpFunctor sort_op)
  {
    cctbx_assert(miller_indices.size() == data.size());
    typedef indexed_value<Index, DataType, SortCmpFunctor> ivalue_type;
    af::shared<ivalue_type> ivalues;
    ivalues.reserve(miller_indices.size());
    std::size_t i;
    for(i=0;i<miller_indices.size();i++) {
      ivalues.push_back(ivalue_type(miller_indices[i], data[i]));
    }
    std::sort(ivalues.begin(), ivalues.end());
    for(i=0;i<miller_indices.size();i++) {
      miller_indices[i] = ivalues[i].index;
      data[i] = ivalues[i].value;
    }
  }

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_SORT_H
