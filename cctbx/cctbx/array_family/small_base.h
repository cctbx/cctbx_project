// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_SMALL_BASE_H
#define CCTBX_ARRAY_FAMILY_SMALL_BASE_H

#include <cctbx/array_family/ref.h>

namespace cctbx { namespace af {

  // Automatic allocation, variable size.
  template <typename ElementType, std::size_t N>
  class small_base
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      ElementType elems[N];

      small_base() : m_size(0) {}

      explicit small_base(size_type n) : m_size(n) {}

      small_base(size_type n, const ElementType& x) : m_size(n) {
        std::fill(begin(), end(), x);
      }

      // Copy with type conversion.
      template <typename OtherElementType, std::size_t OtherN>
      small_base(const small_base<OtherElementType,OtherN>& a)
        : m_size(a.size()) {
        sizecheck(m_size, N);
        for(size_type i=0;i<size();i++) elems[i] = ElementType(a[i]);
      }

      // Assignment with type conversion.
      template <typename OtherElementType, std::size_t OtherN>
      small_base<ElementType,N>&
      operator=(const small_base<OtherElementType,OtherN>& rhs) {
        resize(rhs.size());
        for(size_type i=0;i<size();i++) elems[i] = ElementType(rhs[i]);
        return *this;
      }

      CCTBX_ARRAY_FAMILY_TAKE_REF(elems, N)

      size_type size() const { return m_size; }

      CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(elems, m_size)

      bool empty() const { if (size() == 0) return true; return false; }
      static size_type max_size() { return N; }
      static size_type capacity() { return N; }

      void resize(size_type new_size) {
        sizecheck(new_size, N);
        m_size = new_size;
      }

      void resize(size_type new_size, const ElementType& x) {
        sizecheck(new_size, N);
        if (new_size > m_size) std::fill(end(), begin()+new_size, x);
        m_size = new_size;
      }

      void clear() { m_size = 0; }

      void push_back(const ElementType& x) { elems[m_size++] = x; }
      void pop_back() { m_size--; }

      ElementType* insert(ElementType* pos, size_type n, const ElementType& x) {
        if (n > 0) {
          for (ElementType* p=end()+n-1; p != pos; p--) *p = *(p-n);
          for(size_type i=0;i<n;i++) pos[i] = x;
          m_size += n;
        }
        return pos;
      }

      ElementType* insert(ElementType* pos, const ElementType& x) {
        return insert(pos, 1, x);
      }

      template <typename OtherElementType>
      ElementType* insert(
        ElementType* pos,
        const OtherElementType* first,
        const OtherElementType* last) {
        size_type n = last - first;
        if (n > 0) {
          for (ElementType* p=end()+n-1; p != pos; p--) *p = *(p-n);
          for(size_type i=0;i<n;i++) pos[i] = ElementType(first[i]);
          m_size += n;
        }
        return pos;
      }

      ElementType* erase(ElementType* first, ElementType* last) {
        size_type n = last - first;
        for(ElementType* p=first;p!=end()-n;p++) p[0] = p[n];
        m_size -= n;
        return first;
      }

      ElementType* erase(ElementType* pos) {
        return erase(pos, pos+1);
      }

      void swap(ElementType* other) {
        std::swap(elems, other);
      }

      void assign(size_type n, const ElementType& x) {
        std::fill(begin(), begin()+n, x);
        m_size = n;
      }

      template <typename OtherElementType>
      void assign(
        const OtherElementType* first,
        const OtherElementType* last) {
        m_size = last - first;
        for(size_type i=0;i<size();i++) elems[i] = ElementType(first[i]);
      }

    protected:
      size_type m_size;
  };

}} // namespace cctbx::array_family

#endif // CCTBX_ARRAY_FAMILY_SMALL_BASE_H
