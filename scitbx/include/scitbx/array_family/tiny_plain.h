/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Jan: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_TINY_PLAIN_H
#define SCITBX_ARRAY_FAMILY_TINY_PLAIN_H

#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/detail/misc.h>
#include <scitbx/array_family/detail/tiny_helpers.h>
#include <scitbx/array_family/array_adaptor.h>

namespace scitbx { namespace af {

  // Automatic allocation, fixed size.
  template <typename ElementType, std::size_t N>
  class tiny_plain
  {
    public:
      SCITBX_ARRAY_FAMILY_TYPEDEFS

      ElementType elems[N];

      tiny_plain() {}

      template <typename OtherArrayType>
      tiny_plain(array_adaptor<OtherArrayType> const& a_a)
      {
        OtherArrayType const& a = *(a_a.pointee);
        if (a.size() != N) throw_range_error();
        for(std::size_t i=0;i<N;i++) elems[i] = a[i];
      }

      SCITBX_ARRAY_FAMILY_TINY_CONVENIENCE_CONSTRUCTORS(tiny_plain)
      SCITBX_ARRAY_FAMILY_TINY_COPY_AND_ASSIGNMENT(tiny_plain)

      static size_type size() { return N; }
      static bool empty() { return false; }
      static size_type max_size() { return N; }
      static size_type capacity() { return N; }

      SCITBX_ARRAY_FAMILY_BEGIN_END_ETC(tiny_plain, elems, N)

      SCITBX_ARRAY_FAMILY_TAKE_REF(elems, N)

      void swap(ElementType* other) {
        std::swap(*this, other);
      }
  };

  template <typename ArrayType>
  const_ref<typename ArrayType::value_type>
  make_const_ref(ArrayType const& a) {
    typedef typename ArrayType::value_type value_type;
    typedef const_ref<value_type> return_type;
    typedef typename return_type::accessor_type accessor_type;
    return return_type(&(*(a.begin())), accessor_type(a.size()));
  }

  template <typename ArrayType>
  ref<typename ArrayType::value_type>
  make_ref(ArrayType& a) {
    typedef typename ArrayType::value_type value_type;
    typedef ref<value_type> return_type;
    typedef typename return_type::accessor_type accessor_type;
    return return_type(&(*(a.begin())), accessor_type(a.size()));
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_TINY_PLAIN_H
