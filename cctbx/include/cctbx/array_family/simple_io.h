// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_SIMPLE_IO_H
#define CCTBX_ARRAY_FAMILY_SIMPLE_IO_H

#include <iostream>
#include <cctbx/array_family/tiny.h>

namespace cctbx { namespace af {

  template <typename ElementType>
  std::ostream&
  operator<<(std::ostream& os, const const_ref<ElementType>& a) {
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

  template <typename ElementType>
  std::ostream&
  operator<<(std::ostream& os, const ref<ElementType>& a) {
    return os << a.const_ref();
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SIMPLE_IO_H
