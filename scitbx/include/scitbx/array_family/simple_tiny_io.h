/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Mar: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_SIMPLE_TINY_IO_H
#define SCITBX_ARRAY_FAMILY_SIMPLE_TINY_IO_H

#include <scitbx/array_family/simple_io.h>

namespace scitbx { namespace af {

  template <typename ElementType, std::size_t N>
  std::ostream&
  operator<<(std::ostream& os, tiny<ElementType, N> const& a) {
    return os << a.const_ref();
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_SIMPLE_TINY_IO_H
