// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Dec 2001: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_VECREF_H
#define CCTBX_VECREF_H

#include <cstddef>

namespace cctbx {

  template <typename ValueType>
  class vecref
  {
    public:
      typedef ValueType        value_type;
      typedef ValueType*       iterator;
      typedef const ValueType* const_iterator;
      typedef ValueType&       reference;
      typedef const ValueType& const_reference;
      typedef std::size_t      size_type;
      typedef std::ptrdiff_t   difference_type;

      vecref() : m_begin(0), m_size(0) {}
      template <typename VectorType>
      vecref(VectorType& vec)
        : m_begin(&(*(vec.begin()))), m_size(vec.size()) {}
      template <typename IteratorOrPointerType>
      vecref(IteratorOrPointerType begin, size_type size)
        : m_begin(&(*(begin))), m_size(size) {}
      vecref(void* begin, size_type size)
        : m_begin(reinterpret_cast<ValueType*>(begin)), m_size(size) {}

      ValueType* begin() const { return m_begin; }
      ValueType* end() const { return m_begin + m_size; }

      ValueType& operator[](size_type i) const { return m_begin[i]; }

      size_type size() const { return m_size; }

      void* handle() const { return m_begin; }

    protected:
      ValueType* m_begin;
      std::size_t m_size;
  };

}

#endif // CCTBX_VECREF_H
