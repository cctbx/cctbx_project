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
      shared_plain(const size_type& sz = 0)
        : shared_base<ElementType>(sz)
      {}

      explicit shared_plain(const handle_type& handle)
        : shared_base<ElementType>(handle)
      {}

      CCTBX_ARRAY_FAMILY_TAKE_REF(this->begin(), this->size())

      void auto_resize(const size_type& sz) {
        this->handle().auto_resize(this->element_size() * sz);
      }

      void push_back(const ElementType& value) {
        this->auto_resize(this->size()+1);
        this->operator[](this->size()-1)=value;
      }

      void swap(shared_base<ElementType>& other) {
        this->handle().swap(other.handle());
      }

      void assign(size_type n, const ElementType& x = ElementType()) {
        this->resize(n);
        std::fill(this->begin(), this->end(), x);
      }
  };

}} //namespace cctbx::af

#endif //CCTBX_ARRAY_FAMILY_SHARED_PLAIN_H
