/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Feb: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_VERSA_H
#define SCITBX_ARRAY_FAMILY_VERSA_H

#include <scitbx/array_family/versa_plain.h>
#include <scitbx/array_family/ref_reductions.h>

namespace scitbx { namespace af {

  template <typename ElementType,
            typename AccessorType = trivial_accessor,
            typename BaseArrayType = shared_plain<ElementType> >
  class versa : public versa_plain<ElementType, AccessorType, BaseArrayType>
  {
    public:
      typedef versa<ElementType, AccessorType, BaseArrayType> this_type;

      SCITBX_ARRAY_FAMILY_TYPEDEFS

      typedef BaseArrayType base_array_type;
      typedef versa_plain<ElementType, AccessorType, BaseArrayType> base_class;

      typedef AccessorType accessor_type;
      typedef typename accessor_type::index_type index_type;
      typedef typename accessor_type::index_value_type index_value_type;
      typedef versa<ElementType> one_dim_type;
      typedef typename one_dim_type::accessor_type one_dim_accessor_type;

      versa()
      {}

      explicit
      versa(AccessorType const& ac)
        : base_class(ac)
      {}

      versa(AccessorType const& ac, reserve_flag)
        : base_class(ac, reserve_flag())
      {}

      explicit
      versa(index_value_type const& n0)
        : base_class(n0)
      {}

      versa(AccessorType const& ac, ElementType const& x)
        : base_class(ac, x)
      {}

      versa(index_value_type const& n0, ElementType const& x)
        : base_class(n0, x)
      {}

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
      // non-std
      template <typename FunctorType>
      versa(AccessorType const& ac, init_functor<FunctorType> const& ftor)
        : base_class(ac, ftor)
      {}

      // non-std
      template <typename FunctorType>
      versa(index_value_type const& n0, init_functor<FunctorType> const& ftor)
        : base_class(n0, ftor)
      {}
#endif

      versa(base_class const& other)
        : base_class(other)
      {}

      versa(base_class const& other, weak_ref_flag)
        : base_class(other, weak_ref_flag())
      {}

      versa(base_array_type const& other,
            AccessorType const& ac)
        : base_class(other, ac)
      {}

      versa(base_array_type const& other,
            index_value_type const& n0)
        : base_class(other, n0)
      {}

      versa(base_array_type const& other,
            AccessorType const& ac,
            ElementType const& x)
        : base_class(other, ac, x)
      {}

      versa(base_array_type const& other,
            index_value_type const& n0,
            ElementType const& x)
        : base_class(other, n0, x)
      {}

      versa(sharing_handle* other_handle, AccessorType const& ac)
        : base_class(other_handle, ac)
      {}

      versa(sharing_handle* other_handle, index_value_type const& n0)
        : base_class(other_handle, n0)
      {}

      versa(sharing_handle* other_handle, AccessorType const& ac,
            ElementType const& x)
        : base_class(other_handle, ac)
      {}

      versa(sharing_handle* other_handle, index_value_type const& n0,
            ElementType const& x)
        : base_class(other_handle, n0)
      {}

      template <typename OtherArrayType>
      versa(array_adaptor<OtherArrayType> const& a_a)
        : base_class(a_a)
      {}

      one_dim_type as_1d() {
        return one_dim_type(*this, one_dim_accessor_type(this->size()));
      }

      this_type
      deep_copy() const {
        base_array_type c(this->begin(), this->end());
        return this_type(c, this->m_accessor);
      }

      this_type
      weak_ref() const {
        return this_type(*this, weak_ref_flag());
      }

#     include <scitbx/array_family/detail/reducing_boolean_mem_fun.h>
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_VERSA_H
