// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_SHARED_H
#define CCTBX_ARRAY_FAMILY_SHARED_H

#include <cctbx/array_family/shared_base.h>

namespace cctbx { namespace af {

  // Dynamic allocation, shared (data and size), standard operators.
  template <typename ElementType>
  class shared : public shared_base<ElementType>
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      typedef detail::char_block handle_type;

      explicit shared(const size_type& sz = 0)
        : shared_base<ElementType>(sz)
      {}

      explicit shared(const handle_type& handle)
        : shared_base<ElementType>(handle)
      {}

      CCTBX_ARRAY_FAMILY_TAKE_REF(this->begin(), this->size())
  };

}} // namespace cctbx::af

#include <cctbx/array_family/shared_operators.h>

#endif // CCTBX_ARRAY_FAMILY_SHARED_H
