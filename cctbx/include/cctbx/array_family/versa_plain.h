// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_VERSA_PLAIN_H
#define CCTBX_ARRAY_FAMILY_VERSA_PLAIN_H

#include <cctbx/array_family/shared_plain.h>

namespace cctbx { namespace af {

  template <typename ElementType,
            typename AccessorType = grid<1>,
            typename BaseArrayType = shared_plain<ElementType> >
  class versa_plain : public BaseArrayType
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      typedef BaseArrayType base_array_type;
      typedef BaseArrayType base_class;
      typedef typename base_class::handle_type handle_type;

      typedef AccessorType accessor_type;
      typedef typename accessor_type::index_type index_type;
      typedef versa_plain<ElementType> one_dim_type;
      typedef typename one_dim_type::accessor_type one_dim_accessor_type;

      versa_plain()
      {}

      explicit
      versa_plain(const AccessorType& ac)
        : base_class(ac.size1d()),
          m_accessor(ac)
      {}

      versa_plain(const AccessorType& ac, reserve_flag)
        : base_class(ac.size1d(), reserve_flag()),
          m_accessor(ac)
      {}

      explicit
      versa_plain(long n0)
        : base_class(AccessorType(n0).size1d()),
          m_accessor(n0)
      {}

      versa_plain(const AccessorType& ac, const ElementType& x)
        : base_class(ac.size1d(), x),
          m_accessor(ac)
      {}

      versa_plain(long n0, const ElementType& x)
        : base_class(AccessorType(n0).size1d(), x),
          m_accessor(n0)
      {}

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
      // non-std
      template <typename InitFunctorType>
      versa_plain(const AccessorType& ac, InitFunctorType ftor)
        : base_class(ac.size1d(), ftor),
          m_accessor(ac)
      {}

      // non-std
      template <typename InitFunctorType>
      versa_plain(long n0, InitFunctorType ftor)
        : base_class(AccessorType(n0).size1d(), ftor),
          m_accessor(n0)
      {}
#endif

      versa_plain(const versa_plain<ElementType, AccessorType>& other,
                  weak_ref_flag)
        : base_class(other, weak_ref_flag()),
          m_accessor(other.m_accessor)
      {}

      versa_plain(const base_class& other,
                  const AccessorType& ac)
        : base_class(other),
          m_accessor(ac)
      {
        if (other.size() < size()) throw_range_error();
      }

      versa_plain(const base_class& other,
                  long n0)
        : base_class(other),
          m_accessor(n0)
      {
        if (other.size() < size()) throw_range_error();
      }

      versa_plain(const base_class& other,
                  const AccessorType& ac,
                  const ElementType& x)
        : base_class(other),
          m_accessor(ac)
      {
        base_class::resize(m_accessor.size1d(), x);
      }

      versa_plain(const base_class& other,
                  long n0,
                  const ElementType& x)
        : base_class(other),
          m_accessor(n0)
      {
        base_class::resize(m_accessor.size1d(), x);
      }

      versa_plain(handle_type* other_handle, const AccessorType& ac)
        : base_class(other_handle),
          m_accessor(ac)
      {
        if (other_handle->size / this->element_size() < size()) {
          throw_range_error();
        }
      }

      versa_plain(handle_type* other_handle, long n0)
        : base_class(other_handle),
          m_accessor(n0)
      {
        if (other_handle->size / this->element_size() < size()) {
          throw_range_error();
        }
      }

      versa_plain(handle_type* other_handle, const AccessorType& ac,
                  const ElementType& x)
        : base_class(other_handle),
          m_accessor(ac)
      {
        base_class::resize(m_accessor.size1d(), x);
      }

      versa_plain(handle_type* other_handle, long n0,
                  const ElementType& x)
        : base_class(other_handle),
          m_accessor(n0)
      {
        base_class::resize(m_accessor.size1d(), x);
      }

      const AccessorType& accessor() const { return m_accessor; }
      size_type size() const { return m_accessor.size1d(); }

      // since size() is not a virtual function end() needs to be redefined.
      CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(base_class::begin(), size())

      CCTBX_ARRAY_FAMILY_TAKE_VERSA_REF(begin(), m_accessor)

      void resize(const AccessorType& ac) {
        m_accessor = ac;
        base_class::resize(m_accessor.size1d(), ElementType());
      }

      void resize(const AccessorType& ac, const ElementType& x) {
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

            value_type& operator()(const index_type& i)       {
        return begin()[m_accessor(i)];
      }
      const value_type& operator()(const index_type& i) const {
        return begin()[m_accessor(i)];
      }

      // Convenience operator()

      const value_type& operator()(long i0) const {
        return operator()(index_type(i0));
      }
            value_type& operator()(long i0)       {
        return operator()(index_type(i0));
      }
      const value_type& operator()(long i0,
                                   long i1) const {
        return operator()(index_type(i0, i1));
      }
            value_type& operator()(long i0,
                                   long i1)       {
        return operator()(index_type(i0, i1));
      }
      const value_type& operator()(long i0,
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

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_VERSA_PLAIN_H
