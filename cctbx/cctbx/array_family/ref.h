// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_REF_H
#define CCTBX_ARRAY_FAMILY_REF_H

#include <cstddef>
#include <stdexcept>

// FIXES for broken compilers
#include <boost/config.hpp>

// XXX
#include <iostream>
#define CheckPoint std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl << std::flush

#include <cctbx/array_family/ref_helpers.h>

namespace cctbx { namespace af {

  void throw_range_error() {
    throw std::range_error("array_family");
  }

  template <typename InputElementType,
            typename OutputElementType>
  OutputElementType*
  copy_typeconv(
    const InputElementType* first,
    const InputElementType* last,
    OutputElementType* result)
  {
    OutputElementType* p = result;
    while (first != last) *p++ = OutputElementType(*first++);
    return result;
  }

  template <typename ElementType>
  class const_ref
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      const_ref()
        : m_begin(0), m_size(0)
      {}
      const_ref(const ElementType* begin, size_type size)
        : m_begin(begin), m_size(size)
      {}
      const_ref(const void* begin, size_type size)
        : m_begin(reinterpret_cast<const ElementType*>(begin)), m_size(size)
      {}

      size_type size() const { return m_size; }

      const ElementType* begin() const { return m_begin; }
      const ElementType* end() const { return m_begin + m_size; }
      const ElementType& front() const { return m_begin[0]; }
      const ElementType& back() const { return m_begin[m_size-1]; }

      const ElementType& operator[](size_type i) const { return m_begin[i]; }

      const ElementType& at(size_type i) const {
        if (i >= m_size) throw_range_error();
        return m_begin[i];
      }

      const void* memory_handle() const {
        return reinterpret_cast<const void*>(m_begin);
      }

    private:
      const ElementType* m_begin;
      size_type m_size;
  };

  template <typename ElementType>
  class ref
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      ref()
        : m_begin(0), m_size(0)
      {}
      ref(ElementType* begin, size_type size)
        : m_begin(begin), m_size(size)
      {}
      ref(void* begin, size_type size)
        : m_begin(reinterpret_cast<ElementType*>(begin)), m_size(size)
      {}

      size_type size() const { return m_size; }

      ElementType* begin() const { return m_begin; }
      ElementType* end() const { return m_begin + m_size; }
      ElementType& front() const { return m_begin[0]; }
      ElementType& back() const { return m_begin[m_size-1]; }

      ElementType& operator[](size_type i) const { return m_begin[i]; }

      ElementType& at(size_type i) const {
        if (i >= m_size) throw_range_error();
        return m_begin[i];
      }

      void* memory_handle() const {
        return reinterpret_cast<void*>(m_begin);
      }

      operator const_ref<ElementType>() const {
        return const_ref<ElementType>(m_begin, m_size);
      }

    private:
      ElementType* m_begin;
      size_type m_size;
  };

}} // namespace cctbx::array_family

#endif // CCTBX_ARRAY_FAMILY_REF_H
