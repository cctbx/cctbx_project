// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_TINY_PLAIN_H
#define CCTBX_ARRAY_FAMILY_TINY_PLAIN_H

#include <algorithm>
#include <cctbx/array_family/ref.h>
#include <cctbx/array_family/tiny_helpers.h>

namespace cctbx { namespace af {

  // Automatic allocation, fixed size.
  template <typename ElementType, std::size_t N>
  class tiny_plain
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      ElementType elems[N];

      tiny_plain() {}

      CCTBX_ARRAY_FAMILY_TINY_CONVENIENCE_CONSTRUCTORS(tiny_plain)
      CCTBX_ARRAY_FAMILY_TINY_COPY_AND_ASSIGNMENT(tiny_plain)

      static size_type size() { return N; }
      static bool empty() { return false; }
      static size_type max_size() { return N; }
      static size_type capacity() { return N; }

      CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(elems, N)

      CCTBX_ARRAY_FAMILY_TAKE_REF(elems, N)

      void swap(ElementType* other) {
        std::swap(*this, other);
      }
  };

}} // namespace cctbx::af

#include <cctbx/array_family/tiny_plain_apply.h>

#endif // CCTBX_ARRAY_FAMILY_TINY_PLAIN_H
