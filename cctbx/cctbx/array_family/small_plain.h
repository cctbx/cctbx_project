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
#include <cctbx/array_family/auto_allocator.h>

namespace cctbx { namespace af {

  // Automatic allocation, variable size.
  template <typename ElementType, std::size_t N>
  class small_plain
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      CCTBX_ARRAY_FAMILY_SMALL_CONSTRUCTORS(small_plain)
      CCTBX_ARRAY_FAMILY_SMALL_COPY_AND_ASSIGNMENT(small_plain)

      ~small_plain() { clear(); }

      CCTBX_ARRAY_FAMILY_TAKE_REF(begin(), N)

      size_type size() const { return m_size; }
      bool empty() const { if (size() == 0) return true; return false; }
      static size_type max_size() { return N; }
      static size_type capacity() { return N; }

      CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(
        ((ElementType*)(m_elems.buffer)), m_size) // fix this

      void swap(small_plain<ElementType, N>& other) {
        std::swap(*this, other);
      }

      void reserve(const size_type& sz) {
        if (N < sz) throw_range_error();
      }

#     include <cctbx/array_family/push_back_etc.h>

    protected:
      void m_insert_overflow(ElementType* pos,
                             const size_type& n, const ElementType& x,
                             bool at_end) {
        throw_range_error();
      }

      void m_insert_overflow(ElementType* pos,
                             const ElementType* first,
                             const ElementType* last) {
        throw_range_error();
      }

      void m_set_size(const size_type& sz) {
        m_size = sz;
      }

      void m_incr_size(const size_type& n) {
        m_size += n;
      }

      void m_decr_size(const size_type& n) {
        m_size -= n;
      }

      detail::auto_allocator<ElementType, N> m_elems;
      size_type m_size;
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SMALL_PLAIN_H
