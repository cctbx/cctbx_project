/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Feb: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_VERSA_PLAIN_H
#define SCITBX_ARRAY_FAMILY_VERSA_PLAIN_H

#include <scitbx/array_family/shared_plain.h>

namespace scitbx { namespace af {

  template <typename ElementType,
            typename AccessorType = grid<1>,
            typename BaseArrayType = shared_plain<ElementType> >
  class versa_plain : public BaseArrayType
  {
    public:
      SCITBX_ARRAY_FAMILY_TYPEDEFS

      typedef BaseArrayType base_array_type;
      typedef BaseArrayType base_class;

      typedef AccessorType accessor_type;
      typedef typename accessor_type::index_type index_type;
      typedef versa_plain<ElementType> one_dim_type;
      typedef typename one_dim_type::accessor_type one_dim_accessor_type;

      versa_plain()
      {}

      explicit
      versa_plain(AccessorType const& ac)
        : base_class(ac.size1d()),
          m_accessor(ac)
      {}

      versa_plain(AccessorType const& ac, reserve_flag)
        : base_class(ac.size1d(), reserve_flag()),
          m_accessor(ac)
      {}

      explicit
      versa_plain(long n0)
        : base_class(AccessorType(n0).size1d()),
          m_accessor(n0)
      {}

      versa_plain(AccessorType const& ac, ElementType const& x)
        : base_class(ac.size1d(), x),
          m_accessor(ac)
      {}

      versa_plain(long n0, ElementType const& x)
        : base_class(AccessorType(n0).size1d(), x),
          m_accessor(n0)
      {}

      // non-std
      template <typename FunctorType>
      versa_plain(AccessorType const& ac, init_functor<FunctorType> const& ftor)
        : base_class(ac.size1d(), ftor),
          m_accessor(ac)
      {}

      // non-std
      template <typename FunctorType>
      versa_plain(long n0, init_functor<FunctorType> const& ftor)
        : base_class(AccessorType(n0).size1d(), ftor),
          m_accessor(n0)
      {}

      versa_plain(versa_plain<ElementType, AccessorType> const& other,
                  weak_ref_flag)
        : base_class(other, weak_ref_flag()),
          m_accessor(other.m_accessor)
      {}

      versa_plain(base_class const& other,
                  AccessorType const& ac)
        : base_class(other),
          m_accessor(ac)
      {
        if (other.size() < size()) throw_range_error();
      }

      versa_plain(base_class const& other,
                  long n0)
        : base_class(other),
          m_accessor(n0)
      {
        if (other.size() < size()) throw_range_error();
      }

      versa_plain(base_class const& other,
                  AccessorType const& ac,
                  ElementType const& x)
        : base_class(other),
          m_accessor(ac)
      {
        base_class::resize(m_accessor.size1d(), x);
      }

      versa_plain(base_class const& other,
                  long n0,
                  ElementType const& x)
        : base_class(other),
          m_accessor(n0)
      {
        base_class::resize(m_accessor.size1d(), x);
      }

      versa_plain(sharing_handle* other_handle, AccessorType const& ac)
        : base_class(other_handle),
          m_accessor(ac)
      {
        if (other_handle->size / this->element_size() < size()) {
          throw_range_error();
        }
      }

      versa_plain(sharing_handle* other_handle, long n0)
        : base_class(other_handle),
          m_accessor(n0)
      {
        if (other_handle->size / this->element_size() < size()) {
          throw_range_error();
        }
      }

      versa_plain(sharing_handle* other_handle, AccessorType const& ac,
                  ElementType const& x)
        : base_class(other_handle),
          m_accessor(ac)
      {
        base_class::resize(m_accessor.size1d(), x);
      }

      versa_plain(sharing_handle* other_handle, long n0,
                  ElementType const& x)
        : base_class(other_handle),
          m_accessor(n0)
      {
        base_class::resize(m_accessor.size1d(), x);
      }

      AccessorType const& accessor() const { return m_accessor; }
      size_type size() const { return m_accessor.size1d(); }

      // since size() is not a virtual function end() needs to be redefined.
      SCITBX_ARRAY_FAMILY_BEGIN_END_ETC(versa_plain,
        base_class::begin(), size())

      SCITBX_ARRAY_FAMILY_TAKE_VERSA_REF(begin(), m_accessor)

      void resize(AccessorType const& ac) {
        m_accessor = ac;
        base_class::resize(m_accessor.size1d(), ElementType());
      }

      void resize(AccessorType const& ac, ElementType const& x) {
        m_accessor = ac;
        base_class::resize(m_accessor.size1d(), x);
      }

      one_dim_type as_1d() {
        return one_dim_type(*this, one_dim_accessor_type(size()));
      }

      versa_plain<ElementType, AccessorType>
      deep_copy() const {
        BaseArrayType c(begin(), end());
        return versa_plain<ElementType, AccessorType>(c, m_accessor);
      }

      BaseArrayType
      as_base_array() const {
        return BaseArrayType(*this);
      }

      versa_plain<ElementType, AccessorType>
      weak_ref() const {
        return versa_plain<ElementType, AccessorType>(*this, weak_ref_flag());
      }

            value_type& operator()(index_type const& i)       {
        return begin()[m_accessor(i)];
      }
      value_type const& operator()(index_type const& i) const {
        return begin()[m_accessor(i)];
      }

      // Convenience operator()

      value_type const& operator()(long i0) const {
        return operator()(index_type(i0));
      }
            value_type& operator()(long i0)       {
        return operator()(index_type(i0));
      }
      value_type const& operator()(long i0,
                                   long i1) const {
        return operator()(index_type(i0, i1));
      }
            value_type& operator()(long i0,
                                   long i1)       {
        return operator()(index_type(i0, i1));
      }
      value_type const& operator()(long i0,
                                   long i1,
                                   long i2) const {
        return operator()(index_type(i0, i1, i2));
      }
            value_type& operator()(long i0,
                                   long i1,
                                   long i2)       {
        return operator()(index_type(i0, i1, i2));
      }

    protected:
      AccessorType m_accessor;

    private:
      // disable modifiers (push_back_etc.h), except for resize()
      void assign();
      void push_back();
      void append();
      void pop_back();
      void insert();
      void erase();
      void clear();
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_VERSA_PLAIN_H
