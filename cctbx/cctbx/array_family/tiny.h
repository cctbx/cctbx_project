// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_TINY_H
#define CCTBX_ARRAY_FAMILY_TINY_H

#include <cctbx/array_family/tiny_plain.h>

namespace cctbx { namespace af {

  // Automatic allocation, fixed size, standard operators.
  template <typename ElementType, std::size_t N>
  class tiny : public tiny_plain<ElementType, N>
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      tiny() {}

      CCTBX_ARRAY_FAMILY_TINY_CONVENIENCE_CONSTRUCTORS(tiny)
      CCTBX_ARRAY_FAMILY_TINY_COPY_AND_ASSIGNMENT(tiny)
      CCTBX_ARRAY_FAMILY_TAKE_REF(this->elems, N)
  };

}} // namespace cctbx::af

#include <cctbx/array_family/tiny_operators.h>

#endif // CCTBX_ARRAY_FAMILY_TINY_H
