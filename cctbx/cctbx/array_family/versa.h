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

namespace cctbx { namespace af {

  template <typename ElementType,
            typename AccessorType = grid_accessor<1> >
  class versa : public versa_plain<ElementType, AccessorType>
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      typedef typename shared_base<ElementType>::handle_type handle_type;

      typedef AccessorType accessor_type;
      typedef typename accessor_type::index_type index_type;

      CCTBX_ARRAY_FAMILY_VERSA_CONSTRUCTORS(versa)

      versa(const handle_type& handle, const accessor_type& ac)
        : versa_plain<ElementType, AccessorType>(handle, ac)
      {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
      }
      versa(const handle_type& handle, const size_type& sz)
        : versa_plain<ElementType, AccessorType>(handle, AccessorType(sz))
      {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
      }

      versa<ElementType> as_1d() {
        return versa<ElementType>(this->handle(), this->size());
      }

      CCTBX_ARRAY_FAMILY_TAKE_VERSA_REF(this->begin(), this->accessor())
  };

}} // namespace cctbx::af

#include <cctbx/array_family/versa_algebra.h>

#endif // CCTBX_ARRAY_FAMILY_VERSA_H
