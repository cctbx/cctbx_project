// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Dec 2001: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_VECTOR_CARRAY_TRAITS_H
#define CCTBX_VECTOR_CARRAY_TRAITS_H

#include <cctbx/carray.h>

namespace cctbx { namespace vector {

template <class T, class U, std::size_t N>
struct algebra_traits<carray<T, N>, carray<U, N> >
{
  typedef typename algebra_traits<T, U>::promotion_type value_type;
  typedef carray<value_type, N> promotion_type;
  typedef carray<bool, N> bool_type;
  typedef qualifier_v_v type_qualifier;
  template <typename ResultType, typename SizeType>
  static void init(ResultType& v, const SizeType& n) {}
};

template <class T, class U, std::size_t N>
struct algebra_traits<carray<T, N>, U>
{
  typedef typename algebra_traits<T, U>::promotion_type value_type;
  typedef carray<value_type, N> promotion_type;
  typedef carray<bool, N> bool_type;
  typedef qualifier_v_s type_qualifier;
  template <typename ResultType, typename SizeType>
  static void init(ResultType& v, const SizeType& n) {}
};

template <class T, class U, std::size_t N>
struct algebra_traits<T, carray<U, N> >
{
  typedef typename algebra_traits<T, U>::promotion_type value_type;
  typedef carray<value_type, N> promotion_type;
  typedef carray<bool, N> bool_type;
  typedef qualifier_s_v type_qualifier;
  template <typename ResultType, typename SizeType>
  static void init(ResultType& v, const SizeType& n) {}
};

}} // namespace cctbx::vector

#endif // CCTBX_VECTOR_CARRAY_TRAITS_H
