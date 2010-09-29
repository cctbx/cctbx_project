#ifndef SCITBX_ARRAY_FAMILY_SHARED_H
#define SCITBX_ARRAY_FAMILY_SHARED_H

#include <scitbx/array_family/shared_plain.h>
#include <scitbx/array_family/ref_reductions.h>

namespace scitbx { namespace af {

  // Dynamic allocation, shared (data and size), standard operators.
  template <typename ElementType>
  class shared : public shared_plain<ElementType>
  {
    public:
      SCITBX_ARRAY_FAMILY_TYPEDEFS

      typedef shared_plain<ElementType> base_class;

      shared()
      {}

      explicit
      shared(size_type const& sz)
        : base_class(sz)
      {}

      // non-std
      shared(af::reserve const& sz)
        : base_class(sz)
      {}

      shared(size_type const& sz, ElementType const& x)
        : base_class(sz, x)
      {}

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
      // non-std
      template <typename FunctorType>
      shared(size_type const& sz, init_functor<FunctorType> const& ftor)
        : base_class(sz, ftor)
      {}
#endif

      shared(const ElementType* first, const ElementType* last)
        : base_class(first, last)
      {}

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
      template <typename OtherElementType>
      shared(const OtherElementType* first, const OtherElementType* last)
        : base_class(first, last)
      {}
#endif

      // non-std
      shared(base_class const& other)
        : base_class(other)
      {}

      // non-std
      shared(base_class const& other, weak_ref_flag)
        : base_class(other, weak_ref_flag())
      {}

      // non-std
      explicit
      shared(sharing_handle* other_handle)
        : base_class(other_handle)
      {}

      // non-std
      shared(sharing_handle* other_handle, weak_ref_flag)
        : base_class(other_handle, weak_ref_flag())
      {}

      // non-std
      template <typename OtherArrayType>
      shared(array_adaptor<OtherArrayType> const& a_a)
        : base_class(a_a)
      {}

      // non-std
      template <class E>
      shared(expression<E> const &e)
        : base_class(e.size(), init_functor_null<ElementType>())
      {
        e.assign_to(this->ref());
      }

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

      /// Expression templates
      //@{
      template <class E>
      shared& operator=(expression<E> const &e) {
        this->ref() = e;
        return *this;
      }

      template <class E>
      shared& operator+=(expression<E> const &e) {
        this->ref() += e;
        return *this;
      }

      template <class E>
      shared& operator-=(expression<E> const &e) {
        this->ref() -= e;
        return *this;
      }

      template <class E>
      shared& operator*=(expression<E> const &e) {
        this->ref() *= e;
        return *this;
      }
      //@}

#     include <scitbx/array_family/detail/reducing_boolean_mem_fun.h>
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_SHARED_H
