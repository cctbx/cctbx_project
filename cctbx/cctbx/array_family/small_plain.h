// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_SMALL_PLAIN_H
#define CCTBX_ARRAY_FAMILY_SMALL_PLAIN_H

#include <algorithm>
#include <cctbx/array_family/ref.h>
#include <cctbx/array_family/small_helpers.h>

namespace cctbx { namespace af {

  // Automatic allocation, variable size.
  template <typename ElementType, std::size_t N>
  class small_plain
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      ElementType elems[N];

      CCTBX_ARRAY_FAMILY_SMALL_CONSTRUCTORS(small_plain)
      CCTBX_ARRAY_FAMILY_SMALL_COPY_AND_ASSIGNMENT(small_plain)
      CCTBX_ARRAY_FAMILY_TAKE_REF(elems, N)

      size_type size() const { return m_size; }

      CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(elems, m_size)

      bool empty() const { if (size() == 0) return true; return false; }
      static size_type max_size() { return N; }
      static size_type capacity() { return N; }

      void resize(size_type new_size) {
        if (new_size > N) throw_range_error();
        m_size = new_size;
      }

      void resize(size_type new_size, const ElementType& x) {
        if (new_size > N) throw_range_error();
        if (new_size > m_size) std::fill(end(), begin()+new_size, x);
        m_size = new_size;
      }

      void clear() { m_size = 0; }

      void push_back(const ElementType& x) { elems[m_size++] = x; }
      void pop_back() { m_size--; }

      ElementType*
      insert(ElementType* pos, size_type n, const ElementType& x) {
        std::copy_backward(pos, end(), end()+n);
        std::fill(pos, pos+n, x);
        m_size += n;
        return pos;
      }

      ElementType*
      insert(ElementType* pos, const ElementType& x) {
        return insert(pos, 1, x);
      }

      template <typename OtherElementType>
      ElementType*
      insert(
        ElementType* pos,
        const OtherElementType* first,
        const OtherElementType* last) {
        size_type n = last - first;
        if (n > 0) {
          std::copy_backward(pos, end(), end()+n);
          copy_typeconv(first, last, pos);
          m_size += n;
        }
        return pos;
      }

      ElementType* erase(ElementType* first, ElementType* last) {
        size_type n = last - first;
        std::copy(last, end(), first);
        m_size -= n;
        return first;
      }

      ElementType* erase(ElementType* pos) {
        return erase(pos, pos+1);
      }

      void swap(ElementType* other) {
        std::swap(elems, other);
      }

      void assign(size_type n, const ElementType& x = ElementType()) {
        std::fill(begin(), begin()+n, x);
        m_size = n;
      }

      template <typename OtherElementType>
      void assign(
        const OtherElementType* first,
        const OtherElementType* last) {
        m_size = last - first;
        copy_typeconv(first, last, elems);
      }

    protected:
      size_type m_size;
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SMALL_PLAIN_H
