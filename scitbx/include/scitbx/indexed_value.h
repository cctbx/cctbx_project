/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Jan: Moved from cctbx to scitbx (rwgk)
     2002 Jun: Created, fragment from cctbx/maps/peak_search.h (rwgk)
 */

#ifndef SCITBX_INDEXED_VALUE_H
#define SCITBX_INDEXED_VALUE_H

#include <functional>

namespace scitbx {

  template <typename IndexType,
            typename ValueType,
            typename SortCmpFunctor = std::less<ValueType> >
  struct indexed_value
  {
    typedef IndexType index_type;
    typedef ValueType value_type;

    indexed_value() {}
    indexed_value(index_type const& i, value_type const& v)
      : index(i), value(v)
    {}

    bool
    operator<(
      indexed_value<IndexType, ValueType, SortCmpFunctor> const& rhs) const
    {
      return SortCmpFunctor()(this->value, rhs.value);
    }

    index_type index;
    value_type value;
  };

} // namespace scitbx

#endif // SCITBX_INDEXED_VALUE_H
