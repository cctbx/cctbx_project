// $Id$
/* Copyright (c) 2002 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (N.K. Sauter)
 */

#ifndef CCTBX_ARRAY_FAMILY_SHARED_PLAIN_H
#define CCTBX_ARRAY_FAMILY_SHARED_PLAIN_H

#include <cctbx/array_family/shared_base.h>

namespace cctbx { namespace af {

  template <typename ElementType>
  class shared_plain : public shared_base<ElementType>
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      typedef typename shared_base<ElementType>::handle_type handle_type;

      explicit
      shared_plain(const size_type& sz = 0,
                   const ElementType& x = ElementType())
        : shared_base<ElementType>(sz, x)
      {}

      template <typename OtherElementType>
      shared_plain(const OtherElementType* first, const OtherElementType* last)
        : shared_base<ElementType>(first, last)
      {}

      explicit shared_plain(const handle_type& handle)
        : shared_base<ElementType>(handle)
      {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
      }

      CCTBX_ARRAY_FAMILY_TAKE_REF(this->begin(), this->size())

      void swap(shared_base<ElementType>& other) {
        this->handle().swap(other.handle());
      }

      void auto_resize(const size_type& sz) {
        this->handle().auto_resize(this->element_size() * sz);
      }

#     include <cctbx/array_family/push_back_etc.h>
  };

}} //namespace cctbx::af

#endif //CCTBX_ARRAY_FAMILY_SHARED_PLAIN_H
