// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_VERSA_H
#define CCTBX_ARRAY_FAMILY_VERSA_H

#include <cctbx/array_family/versa_plain.h>
#include <cctbx/array_family/shared.h>

namespace cctbx { namespace af {

  template <typename ElementType,
            typename AccessorType = grid<1> >
  class versa : public versa_plain<ElementType, AccessorType>
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      typedef AccessorType accessor_type;

      typedef versa_plain<ElementType, AccessorType> base_class;
      typedef typename base_class::handle_type handle_type;
      typedef versa<ElementType> one_dim_type;
      typedef typename one_dim_type::accessor_type one_dim_accessor_type;

      versa()
      {}

      explicit
      versa(const AccessorType& ac)
        : base_class(ac)
      {}

      versa(const AccessorType& ac, reserve_flag)
        : base_class(ac, reserve_flag())
      {}

      explicit
      versa(long n0)
        : base_class(n0)
      {}

      versa(const AccessorType& ac, no_initialization_flag)
        : base_class(ac, no_initialization_flag())
      {}

      versa(long n0, no_initialization_flag)
        : base_class(n0, no_initialization_flag())
      {}

      versa(const AccessorType& ac, const ElementType& x)
        : base_class(ac, x)
      {}

      versa(long n0, const ElementType& x)
        : base_class(n0, x)
      {}

      versa(const versa<ElementType, AccessorType>& other, weak_ref_flag)
        : base_class(other, weak_ref_flag())
      {}

      template <typename OtherAccessorType>
      versa(versa<ElementType, OtherAccessorType>& other,
            const AccessorType& ac)
        : base_class(other, ac)
      {}

      template <typename OtherAccessorType>
      versa(versa<ElementType, OtherAccessorType>& other,
                  long n0)
        : base_class(other, n0)
      {}

      template <typename OtherAccessorType>
      versa(versa<ElementType, OtherAccessorType>& other,
                  const AccessorType& ac,
                  const ElementType& x)
        : base_class(other, ac, x)
      {}

      template <typename OtherAccessorType>
      versa(versa<ElementType, OtherAccessorType>& other,
                  long n0,
                  const ElementType& x)
        : base_class(other, n0, x)
      {}

      versa(handle_type* other_handle, const AccessorType& ac)
        : base_class(other_handle, ac)
      {}

      versa(handle_type* other_handle, long n0)
        : base_class(other_handle, n0)
      {}

      versa(const handle_type& other_handle, const AccessorType& ac,
                  const ElementType& x)
        : base_class(other_handle, ac)
      {}

      versa(const handle_type& other_handle, long n0,
                  const ElementType& x)
        : base_class(other_handle, n0)
      {}

      one_dim_type as_1d() {
        return one_dim_type(*this, one_dim_accessor_type(this->size()));
      }

      versa<ElementType, AccessorType>
      deep_copy() const {
        shared_plain<ElementType> c(this->begin(), this->end());
        return versa<ElementType, AccessorType>(c.handle(), this->m_accessor);
      }

      shared<ElementType>
      as_shared() const {
        return shared<ElementType>(*this);
      }

      versa<ElementType, AccessorType>
      weak_ref() const {
        return versa<ElementType, AccessorType>(*this, weak_ref_flag());
      }
  };

}} // namespace cctbx::af

#include <cctbx/array_family/versa_apply.h>

#endif // CCTBX_ARRAY_FAMILY_VERSA_H
