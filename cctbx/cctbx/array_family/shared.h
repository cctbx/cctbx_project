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

      typedef shared_plain<ElementType> base_class;
      typedef typename base_class::handle_type handle_type;

      shared()
      {}

      explicit
      shared(const size_type& sz)
        : base_class(sz)
      {}

      // non-std
      shared(const size_type& sz, reserve_flag)
        : base_class(sz, reserve_flag())
      {}

      // non-std
      shared(const size_type& sz, no_initialization_flag)
        : base_class(sz, no_initialization_flag())
      {}

      shared(const size_type& sz, const ElementType& x)
        : base_class(sz, x)
      {}

      shared(const ElementType* first, const ElementType* last)
        : base_class(first, last)
      {}

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
      template <typename OtherElementType>
      shared(const OtherElementType* first, const OtherElementType* last)
        : base_class(first, last)
      {}
#endif

      // non-std: hook for versa::as_shared()
      shared(const shared_plain<ElementType>& other)
        : base_class(other)
      {}

      // non-std: shallow copy semantics, weak reference
      shared(const shared<ElementType>& other, weak_ref_flag)
        : base_class(other, weak_ref_flag())
      {}

      // non-std
      explicit
      shared(const handle_type& other_handle)
        : base_class(other_handle)
      {}

      // non-std
      explicit
      shared(const handle_type& other_handle, weak_ref_flag)
        : base_class(other_handle, weak_ref_flag())
      {}

      // non-std
      shared<ElementType>
      deep_copy() const {
        return shared<ElementType>(this->begin(), this->end());
      }

      // non-std
      shared<ElementType>
      weak_ref() const {
        return shared<ElementType>(*this, weak_ref_flag());
      }
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SHARED_H
