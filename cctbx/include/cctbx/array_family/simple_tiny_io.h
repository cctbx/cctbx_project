// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_SIMPLE_TINY_IO_H
#define CCTBX_ARRAY_FAMILY_SIMPLE_TINY_IO_H

#include <cctbx/array_family/simple_io.h>

namespace cctbx { namespace af {

  template <typename ElementType, std::size_t N>
  std::ostream&
  operator<<(std::ostream& os, tiny<ElementType, N> const& a) {
    return os << a.const_ref();
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SIMPLE_TINY_IO_H
