/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Feb: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_SIMPLE_IO_H
#define SCITBX_ARRAY_FAMILY_SIMPLE_IO_H

#include <iostream>
#include <scitbx/array_family/tiny.h>

namespace scitbx { namespace af {

  template <typename ElementType>
  std::ostream&
  operator<<(std::ostream& os, const_ref<ElementType> const& a) {
    os << "(";
    if (a.size() > 0) {
      for (std::size_t i = 0;;) {
        os << a[i];
        i++;
        if (i == a.size()) break;
        os << ",";
      }
    }
    os << ")";
    return os;
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_SIMPLE_IO_H
