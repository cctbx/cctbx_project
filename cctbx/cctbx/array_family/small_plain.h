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
      bool empty() const { if (size() == 0) return true; return false; }
      static size_type max_size() { return N; }
      static size_type capacity() { return N; }

      CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(elems, m_size)

      void swap(small_plain<ElementType, N>& other) {
        std::swap(*this, other);
      }

      void resize(size_type new_size) {
        if (new_size > N) throw_range_error();
        m_size = new_size;
      }

      void auto_resize(size_type new_size) { resize(new_size); }

      void resize(const size_type& new_size, const ElementType& x) {
        ElementType* old_end = this->end();
        this->resize(new_size);
        if (this->end() > old_end) std::fill(old_end, this->end(), x);
      }

#     include <cctbx/array_family/push_back_etc.h>

    protected:
      size_type m_size;
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SMALL_PLAIN_H
