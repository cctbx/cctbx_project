#ifndef SCITBX_ARRAY_FAMILY_VERSA_PLAIN_H
#define SCITBX_ARRAY_FAMILY_VERSA_PLAIN_H

#include <scitbx/error.h>
#include <scitbx/array_family/shared_plain.h>

namespace scitbx { namespace af {

  template <typename ElementType,
            typename AccessorType = trivial_accessor>
  class versa_plain : public shared_plain<ElementType>
  {
    public:
      SCITBX_ARRAY_FAMILY_TYPEDEFS

      typedef shared_plain<ElementType> base_class;
      typedef base_class base_array_type;

      typedef AccessorType accessor_type;
      typedef typename accessor_type::index_type index_type;
      typedef typename accessor_type::index_value_type index_value_type;
      typedef versa_plain<ElementType> one_dim_type;
      typedef typename one_dim_type::accessor_type one_dim_accessor_type;

      versa_plain()
      {}

      explicit
      versa_plain(AccessorType const& ac)
        : base_class(ac.size_1d()),
          m_accessor(ac)
      {}

      explicit
      versa_plain(index_value_type const& n0)
        : base_class(AccessorType(n0).size_1d()),
          m_accessor(n0)
      {}

      versa_plain(AccessorType const& ac, ElementType const& x)
        : base_class(ac.size_1d(), x),
          m_accessor(ac)
      {}

      versa_plain(index_value_type const& n0, ElementType const& x)
        : base_class(AccessorType(n0).size_1d(), x),
          m_accessor(n0)
      {}

      // non-std
      template <typename FunctorType>
      versa_plain(AccessorType const& ac,
                  init_functor<FunctorType> const& ftor)
        : base_class(ac.size_1d(), ftor),
          m_accessor(ac)
      {}

      // non-std
      template <typename FunctorType>
      versa_plain(index_value_type const& n0,
                  init_functor<FunctorType> const& ftor)
        : base_class(AccessorType(n0).size_1d(), ftor),
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
                  index_value_type const& n0)
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
        base_class::resize(m_accessor.size_1d(), x);
      }

      versa_plain(base_class const& other,
                  index_value_type const& n0,
                  ElementType const& x)
        : base_class(other),
          m_accessor(n0)
      {
        base_class::resize(m_accessor.size_1d(), x);
      }

      versa_plain(sharing_handle* other_handle, AccessorType const& ac)
        : base_class(other_handle),
          m_accessor(ac)
      {
        if (other_handle->size / this->element_size() < size()) {
          throw_range_error();
        }
      }

      versa_plain(sharing_handle* other_handle, index_value_type const& n0)
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
        base_class::resize(m_accessor.size_1d(), x);
      }

      versa_plain(sharing_handle* other_handle, index_value_type const& n0,
                  ElementType const& x)
        : base_class(other_handle),
          m_accessor(n0)
      {
        base_class::resize(m_accessor.size_1d(), x);
      }

      template <typename OtherArrayType>
      versa_plain(array_adaptor<OtherArrayType> const& a_a)
        : base_class(a_a),
          m_accessor((a_a.pointee)->size())
      {}

      AccessorType const& accessor() const { return m_accessor; }

      bool check_shared_size() const
      {
        return base_class::size() >= m_accessor.size_1d();
      }

      size_type size() const
      {
        size_type sz = m_accessor.size_1d();
        SCITBX_ASSERT(base_class::size() >= sz);
        return sz;
      }

      // since size() is not a virtual function end() needs to be redefined.
      SCITBX_ARRAY_FAMILY_BEGIN_END_ETC(versa_plain,
        base_class::begin(), size())

      SCITBX_ARRAY_FAMILY_TAKE_VERSA_REF(begin(), m_accessor)

      void resize(AccessorType const& ac) {
        m_accessor = ac;
        base_class::resize(m_accessor.size_1d(), ElementType());
      }

      void resize(AccessorType const& ac, ElementType const& x) {
        m_accessor = ac;
        base_class::resize(m_accessor.size_1d(), x);
      }

      one_dim_type as_1d() {
        return one_dim_type(*this, one_dim_accessor_type(size()));
      }

      versa_plain<ElementType, AccessorType>
      deep_copy() const {
        base_array_type c(begin(), end());
        return versa_plain<ElementType, AccessorType>(c, m_accessor);
      }

      base_array_type
      as_base_array() const {
        return base_array_type(*this);
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

      value_type const& operator()(index_value_type const& i0) const {
        return begin()[m_accessor(i0)];
      }
            value_type& operator()(index_value_type const& i0)       {
        return begin()[m_accessor(i0)];
      }
      value_type const& operator()(index_value_type const& i0,
                                   index_value_type const& i1) const {
        return begin()[m_accessor(i0, i1)];
      }
            value_type& operator()(index_value_type const& i0,
                                   index_value_type const& i1)       {
        return begin()[m_accessor(i0, i1)];
      }
      value_type const& operator()(index_value_type const& i0,
                                   index_value_type const& i1,
                                   index_value_type const& i2) const {
        return begin()[m_accessor(i0, i1, i2)];
      }
            value_type& operator()(index_value_type const& i0,
                                   index_value_type const& i1,
                                   index_value_type const& i2)       {
        return begin()[m_accessor(i0, i1, i2)];
      }
      value_type const& operator()(index_value_type const& i0,
                                   index_value_type const& i1,
                                   index_value_type const& i2,
                                   index_value_type const& i3) const {
        return begin()[m_accessor(i0, i1, i2, i3)];
      }
            value_type& operator()(index_value_type const& i0,
                                   index_value_type const& i1,
                                   index_value_type const& i2,
                                   index_value_type const& i3)       {
        return begin()[m_accessor(i0, i1, i2, i3)];
      }
      value_type const& operator()(index_value_type const& i0,
                                   index_value_type const& i1,
                                   index_value_type const& i2,
                                   index_value_type const& i3,
                                   index_value_type const& i4) const {
        return begin()[m_accessor(i0, i1, i2, i3, i4)];
      }
            value_type& operator()(index_value_type const& i0,
                                   index_value_type const& i1,
                                   index_value_type const& i2,
                                   index_value_type const& i3,
                                   index_value_type const& i4)       {
        return begin()[m_accessor(i0, i1, i2, i3, i4)];
      }

    protected:
      AccessorType m_accessor;

    private:
      // disable modifiers (push_back_etc.h), except for resize()
      void assign();
      void push_back();
      void append();
      void extend();
      void pop_back();
      void insert();
      void erase();
      void clear();
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_VERSA_PLAIN_H
