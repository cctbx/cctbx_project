// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SHARED_STORAGE_H
#define CCTBX_SHARED_STORAGE_H

#include <cstddef>
#include <boost/smart_ptr.hpp>

namespace cctbx {

  template <typename ValueType>
  class shared_storage
  {
    public:
      typedef ValueType        value_type;
      typedef ValueType*       iterator;
      typedef const ValueType* const_iterator;
      typedef ValueType&       reference;
      typedef const ValueType& const_reference;
      typedef std::size_t      size_type;
      typedef std::ptrdiff_t   difference_type;

      typedef boost::shared_array<char> shared_array_type;

      explicit shared_storage(const size_type& sz = 0)
        : m_begin(new value_type[sz]),
          m_size(sz)
      {
        char* char_ptr_begin = reinterpret_cast<char*>(m_begin);
        m_storage = shared_array_type(char_ptr_begin);
      }
      shared_storage(const shared_array_type& storage, size_type size)
        : m_storage(storage),
          m_begin(reinterpret_cast<value_type*>(storage.get())),
          m_size(size)
      {}

      value_type* begin() const { return m_begin; }
      value_type* end() const { return m_begin + m_size; }

      value_type& operator[](size_type i) const { return m_begin[i]; }

      size_type size() const { return m_size; }

      shared_array_type cast() const { return m_storage; }

      shared_storage<value_type>
      deepcopy() const
      {
        shared_storage<value_type> result(m_size);
        std::copy(this->begin(), this->end(), result->begin());
        return result;
      }

    protected:
      shared_array_type m_storage;
      value_type* m_begin;
      size_type m_size;
  };

}

#endif // CCTBX_SHARED_STORAGE_H
