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

#include <cctbx/array_family/apply.h>
#include <cctbx/array_family/shared_plain.h>

namespace cctbx { namespace af {

  // Dynamic allocation, shared (data and size), standard operators.
  template <typename ElementType>
  class shared : public shared_plain<ElementType>
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      typedef typename shared_base<ElementType>::handle_type handle_type;

      explicit
      shared(const size_type& sz = 0,
             const ElementType& x = ElementType())
        : shared_plain<ElementType>(sz, x)
      {}

      template <typename OtherElementType>
      shared(const OtherElementType* first, const OtherElementType* last)
        : shared_plain<ElementType>(first, last)
      {}

      explicit shared(const handle_type& handle)
        : shared_plain<ElementType>(handle)
      {}

      CCTBX_ARRAY_FAMILY_TAKE_REF(this->begin(), this->size())
  };

  template <typename ElementType, typename OtherElementType>
  struct change_array_element_type<
    shared<ElementType>, OtherElementType> {
    typedef shared<OtherElementType> array_type;
  };

}} // namespace cctbx::af

#include <cctbx/array_family/shared_algebra.h>

#endif // CCTBX_ARRAY_FAMILY_SHARED_H
