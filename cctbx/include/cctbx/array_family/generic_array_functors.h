// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_GENERIC_ARRAY_FUNCTORS_H
#define CCTBX_ARRAY_FAMILY_GENERIC_ARRAY_FUNCTORS_H

#include <cctbx/array_family/generic_array_operators.h>

namespace cctbx { namespace af {

  template <typename FunctorType,
            typename ElementType,
            typename ElementTypeResult>
  class array_functor_a {
    public:
      array_functor_a(FunctorType const& ftor, const ElementType* a)
        : m_ftor(ftor), m_a(a)
      {}
      void operator()(
        ElementTypeResult* result,
        std::size_t const& sz) const {
        array_operation_a(m_ftor, m_a, result, sz,
          has_trivial_destructor<ElementTypeResult>::value());
      }
    protected:
      FunctorType const& m_ftor;
      const ElementType* m_a;
  };

  template <typename FunctorType,
            typename ElementType>
  inline
  array_functor_a<
    FunctorType,
    ElementType,
    typename FunctorType::result_type>
  make_array_functor_a(
    FunctorType const& ftor,
    const ElementType* a) {
      return array_functor_a<
        FunctorType,
        ElementType,
        typename FunctorType::result_type>(ftor, a);
  }

  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  class array_functor_a_a {
    public:
      array_functor_a_a(
       FunctorType const& ftor,
       const ElementType1* a1,
       const ElementType2* a2)
        : m_ftor(ftor), m_a1(a1), m_a2(a2)
      {}
      void operator()(
        ElementTypeResult* result,
        std::size_t const& sz) const {
        array_operation_a_a(m_ftor, m_a1, m_a2, result, sz,
          has_trivial_destructor<ElementTypeResult>::value());
      }
    protected:
      FunctorType const& m_ftor;
      const ElementType1* m_a1;
      const ElementType2* m_a2;
  };

  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2>
  inline
  array_functor_a_a<
    FunctorType,
    ElementType1,
    ElementType2,
    typename FunctorType::result_type>
  make_array_functor_a_a(
    FunctorType const& ftor,
    const ElementType1* a1,
    const ElementType2* a2) {
      return array_functor_a_a<
        FunctorType,
        ElementType1,
        ElementType2,
        typename FunctorType::result_type>(ftor, a1, a2);
  }

  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  class array_functor_a_s {
    public:
      array_functor_a_s(
       FunctorType const& ftor,
       const ElementType1* a1,
       ElementType2 const& a2)
        : m_ftor(ftor), m_a1(a1), m_a2(a2)
      {}
      void operator()(
        ElementTypeResult* result,
        std::size_t const& sz) const {
        array_operation_a_s(m_ftor, m_a1, m_a2, result, sz,
          has_trivial_destructor<ElementTypeResult>::value());
      }
    protected:
      FunctorType const& m_ftor;
      const ElementType1* m_a1;
      ElementType2 m_a2;
  };

  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2>
  inline
  array_functor_a_s<
    FunctorType,
    ElementType1,
    ElementType2,
    typename FunctorType::result_type>
  make_array_functor_a_s(
    FunctorType const& ftor,
    const ElementType1* a1,
    ElementType2 const& a2) {
      return array_functor_a_s<
        FunctorType,
        ElementType1,
        ElementType2,
        typename FunctorType::result_type>(ftor, a1, a2);
  }

  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  class array_functor_s_a {
    public:
      array_functor_s_a(
       FunctorType const& ftor,
       ElementType1 const& a1,
       const ElementType2* a2)
        : m_ftor(ftor), m_a1(a1), m_a2(a2)
      {}
      void operator()(
        ElementTypeResult* result,
        std::size_t const& sz) const {
        array_operation_s_a(m_ftor, m_a1, m_a2, result, sz,
          has_trivial_destructor<ElementTypeResult>::value());
      }
    protected:
      FunctorType const& m_ftor;
      ElementType1 m_a1;
      const ElementType2* m_a2;
  };

  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2>
  inline
  array_functor_s_a<
    FunctorType,
    ElementType1,
    ElementType2,
    typename FunctorType::result_type>
  make_array_functor_s_a(
    FunctorType const& ftor,
    ElementType1 const& a1,
    const ElementType2* a2) {
      return array_functor_s_a<
        FunctorType,
        ElementType1,
        ElementType2,
        typename FunctorType::result_type>(ftor, a1, a2);
  }

  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementType3,
            typename ElementTypeResult>
  class array_functor_a_a_s {
    public:
      array_functor_a_a_s(
       FunctorType const& ftor,
       const ElementType1* a1,
       const ElementType2* a2,
       ElementType3 const& a3)
        : m_ftor(ftor), m_a1(a1), m_a2(a2), m_a3(a3)
      {}
      void operator()(
        ElementTypeResult* result,
        std::size_t const& sz) const {
        array_operation_a_a_s(m_ftor, m_a1, m_a2, m_a3, result, sz,
          has_trivial_destructor<ElementTypeResult>::value());
      }
    protected:
      FunctorType const& m_ftor;
      const ElementType1* m_a1;
      const ElementType2* m_a2;
      ElementType3 m_a3;
  };

  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementType3>
  inline
  array_functor_a_a_s<
    FunctorType,
    ElementType1,
    ElementType2,
    ElementType3,
    typename FunctorType::result_type>
  make_array_functor_a_a_s(
    FunctorType const& ftor,
    const ElementType1* a1,
    const ElementType2* a2,
    ElementType3 const& a3) {
      return array_functor_a_a_s<
        FunctorType,
        ElementType1,
        ElementType2,
        ElementType3,
        typename FunctorType::result_type>(ftor, a1, a2, a3);
  }

  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementType3,
            typename ElementTypeResult>
  class array_functor_a_s_s {
    public:
      array_functor_a_s_s(
       FunctorType const& ftor,
       const ElementType1* a1,
       ElementType2 const& a2,
       ElementType3 const& a3)
        : m_ftor(ftor), m_a1(a1), m_a2(a2), m_a3(a3)
      {}
      void operator()(
        ElementTypeResult* result,
        std::size_t const& sz) const {
        array_operation_a_s_s(m_ftor, m_a1, m_a2, m_a3, result, sz,
          has_trivial_destructor<ElementTypeResult>::value());
      }
    protected:
      FunctorType const& m_ftor;
      const ElementType1* m_a1;
      ElementType2 m_a2;
      ElementType3 m_a3;
  };

  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementType3>
  inline
  array_functor_a_s_s<
    FunctorType,
    ElementType1,
    ElementType2,
    ElementType3,
    typename FunctorType::result_type>
  make_array_functor_a_s_s(
    FunctorType const& ftor,
    const ElementType1* a1,
    ElementType2 const& a2,
    ElementType3 const& a3) {
      return array_functor_a_s_s<
        FunctorType,
        ElementType1,
        ElementType2,
        ElementType3,
        typename FunctorType::result_type>(ftor, a1, a2, a3);
  }

  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementType3,
            typename ElementTypeResult>
  class array_functor_s_a_s {
    public:
      array_functor_s_a_s(
       FunctorType const& ftor,
       ElementType1 const& a1,
       const ElementType2* a2,
       ElementType3 const& a3)
        : m_ftor(ftor), m_a1(a1), m_a2(a2), m_a3(a3)
      {}
      void operator()(
        ElementTypeResult* result,
        std::size_t const& sz) const {
        array_operation_s_a_s(m_ftor, m_a1, m_a2, m_a3, result, sz,
          has_trivial_destructor<ElementTypeResult>::value());
      }
    protected:
      FunctorType const& m_ftor;
      ElementType1 m_a1;
      const ElementType2* m_a2;
      ElementType3 m_a3;
  };

  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementType3>
  inline
  array_functor_s_a_s<
    FunctorType,
    ElementType1,
    ElementType2,
    ElementType3,
    typename FunctorType::result_type>
  make_array_functor_s_a_s(
    FunctorType const& ftor,
    ElementType1 const& a1,
    const ElementType2* a2,
    ElementType3 const& a3) {
      return array_functor_s_a_s<
        FunctorType,
        ElementType1,
        ElementType2,
        ElementType3,
        typename FunctorType::result_type>(ftor, a1, a2, a3);
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_GENERIC_ARRAY_FUNCTORS_H
