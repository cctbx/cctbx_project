// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_BASIC_ARRAY_FAMILY_H
#define CCTBX_BASIC_ARRAY_FAMILY_H

#include <cstddef>
#include <stdexcept>
#include <algorithm>

// FIXES for broken compilers
#include <boost/config.hpp>

#include <array_family_helpers.h>

namespace cctbx { namespace array_family {

  void rangecheck (size_type i, size_type n) {
    if (i >= n) { throw std::range_error("array_family"); }
  }
  void sizecheck (size_type s, size_type n) {
    if (s > n) { throw std::range_error("array_family"); }
  }

  template <typename ElementType>
  class const_ref
  {
    public:
      CCTBX_BASIC_ARRAY_FAMILY_TYPEDEFS

      const_ref()
        : m_begin(0), m_size(0)
      {}
      const_ref(const ElementType* begin, size_type size)
        : m_begin(begin), m_size(size)
      {}
      const_ref(const void* begin, size_type size)
        : m_begin(reinterpret_cast<const ElementType*>(begin)), m_size(size)
      {}

      const ElementType* begin() const { return m_begin; }
      const ElementType* end() const { return m_begin + m_size; }
      const ElementType& front() const { return m_begin[0]; }
      const ElementType& back() const { return m_begin[size()-1]; }

      const ElementType& operator[](size_type i) const { return m_begin[i]; }
      size_type size() const { return m_size; }

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
      CCTBX_BASIC_ARRAY_FAMILY_TYPEDEFS

      ref()
        : m_begin(0), m_size(0)
      {}
      ref(ElementType* begin, size_type size)
        : m_begin(begin), m_size(size)
      {}
      ref(void* begin, size_type size)
        : m_begin(reinterpret_cast<ElementType*>(begin)), m_size(size)
      {}

      ElementType* begin() const { return m_begin; }
      ElementType* end() const { return m_begin + m_size; }
      ElementType& front() const { return m_begin[0]; }
      ElementType& back() const { return m_begin[size()-1]; }

      ElementType& operator[](size_type i) const { return m_begin[i]; }

      void* memory_handle() const {
        return reinterpret_cast<void*>(m_begin);
      }

      operator const_ref<ElementType>() const {
        const_ref<ElementType>()(m_begin, m_size);
      }

    private:
      ElementType* m_begin;
      size_type m_size;
  };

  // Automatic allocation, fixed size.
  template <typename ElementType, std::size_t N>
  class tiny
  {
    public:
      CCTBX_BASIC_ARRAY_FAMILY_TYPEDEFS

      ElementType elems[N];

      tiny() {}

      // Convenience constructors.
      CCTBX_BASIC_ARRAY_FAMILY_CONVENIENCE_CONSTRUCTORS_FIXSIZE(tiny)

      // Copy with type conversion.
      template <typename OtherElementType>
      tiny(const tiny<OtherElementType,N>& a) {
        for(std::size_t i=0;i<size();i++) elems[i] = a[i];
      }

      // Assignment with type conversion.
      template <typename OtherElementType>
      tiny<ElementType,N>& operator=(const tiny<OtherElementType,N>& rhs) {
        std::copy(rhs.begin(),rhs.end(), begin());
        return *this;
      }

      ElementType* begin() { return elems; }
      const ElementType* begin() const { return elems; }
      ElementType* end() { return elems+N; }
      const ElementType* end() const { return elems+N; }
      ElementType& front() { return elems[0]; }
      const ElementType& front() const { return elems[0]; }
      ElementType& back() { return elems[N-1]; }
      const ElementType& back() const { return elems[N-1]; }

      ElementType& operator[](size_type i) { return elems[i]; }
      const ElementType& operator[](size_type i) const { return elems[i]; }

      ElementType& at(size_type i) {
        rangecheck(i,N); return elems[i];
      }
      const ElementType& at(size_type i) const {
        rangecheck(i,N); return elems[i];
      }

      static size_type size() { return N; }
      static bool empty() { return false; }
      static size_type max_size() { return N; }
      static size_type capacity() { return N; }
  };

  // Automatic allocation, variable size.
  template <typename ElementType, std::size_t N>
  class auto
  {
    public:
      CCTBX_BASIC_ARRAY_FAMILY_TYPEDEFS

      ElementType elems[N];

      auto() : m_size(0) {}

      // Convenience constructors.
      CCTBX_BASIC_ARRAY_FAMILY_CONVENIENCE_CONSTRUCTORS_VARSIZE(auto)

      // Copy with type conversion.
      template <typename OtherElementType, typename OtherN>
      auto(const auto<OtherElementType,OtherN>& a)
        : m_size(a.size()) {
        sizecheck(m_size, N);
        for(std::size_t i=0;i<size();i++) elems[i] = a[i];
      }

      // Assignment with type conversion.
      template <typename OtherElementType, typename OtherN>
      auto<ElementType,N>&
      operator=(const auto<OtherElementType,OtherN>& rhs) {
        sizecheck(rhs.size(), N);
        std::copy(rhs.begin(),rhs.end(), begin());
        return *this;
      }

      ElementType* begin() { return elems; }
      const ElementType* begin() const { return elems; }
      ElementType* end() { return elems+size(); }
      const ElementType* end() const { return elems+size(); }
      ElementType& front() { return elems[0]; }
      const ElementType& front() const { return elems[0]; }
      ElementType& back() { return elems[size()-1]; }
      const ElementType& back() const { return elems[size()-1]; }

      ElementType& operator[](size_type i) { return elems[i]; }
      const ElementType& operator[](size_type i) const { return elems[i]; }

      ElementType& at(size_type i) {
        rangecheck(i,size()); return elems[i];
      }
      const ElementType& at(size_type i) const {
        rangecheck(i,size()); return elems[i];
      }

      size_type size() const { return m_size; }
      bool empty() const { if (size() == 0) return true; return false; }
      static size_type max_size() { return N; }
      static size_type capacity() { return N; }

      void resize (size_type new_size, bool c = false) { // XXX what is c?
        rangecheck(new_size, N);
        m_size = new_size;
      }

      void push_back(const ElementType& x) {
        rangecheck(m_size+1, N);
        elems[m_size] = x;
        m_size++;
      }

/* XXX TODO
      void pop_back () { --_RWfinish; }

      iterator insert (iterator position, const bool& x = bool())
      {
        size_type n = position - begin();
        if (_RWfinish.p != _RWend_of_storage.data() && position == end())
          *_RWfinish++ = x;
        else
          _RWinsert_aux(position, x);
        return begin() + n;
      }

      void insert (iterator position, size_type n, const bool& x);

      template<class InputIterator>
      void insert (iterator position, InputIterator first, InputIterator last);

      iterator erase (iterator position)
      {
        if (!(position + 1 == end()))
          _RWcopy(position + 1, end(), position);
        --_RWfinish;
        return position;
      }

      iterator erase(iterator first, iterator last)
      {
        _RWfinish = _RWcopy(last, end(), first);
        return first;
      }

      void swap (vector<bool,Allocator >& x);

      static void swap(reference x, reference y);

      void flip ();

      void clear()
      {
        erase(begin(),end());
      }
*/

    protected:
      size_type m_size;
  };

}} // namespace cctbx::array_family

#endif CCTBX_BASIC_ARRAY_FAMILY_H
