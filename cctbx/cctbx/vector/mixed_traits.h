// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Dec 2001: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_VECTOR_MIXED_TRAITS_H
#define CCTBX_VECTOR_MIXED_TRAITS_H

#include <cctbx/vector/std_vector_traits.h>
#include <cctbx/vector/carray_traits.h>

namespace cctbx { namespace vector {

template <class T, class U, std::size_t N>
struct algebra_traits<std::vector<T>, carray<U, N> >
{
  typedef typename algebra_traits<T, U>::promotion_type value_type;
  typedef std::vector<value_type> promotion_type;
  typedef std::vector<bool> bool_type;
  typedef qualifier_v_v type_qualifier;
  template <typename ResultType, typename SizeType>
  static void init(ResultType& v, const SizeType& n) { v.resize(n); }
};

template <class T, class U, std::size_t N>
struct algebra_traits<carray<T, N>, std::vector<U> >
{
  typedef typename algebra_traits<T, U>::promotion_type value_type;
  typedef std::vector<value_type> promotion_type;
  typedef std::vector<bool> bool_type;
  typedef qualifier_v_v type_qualifier;
  template <typename ResultType, typename SizeType>
  static void init(ResultType& v, const SizeType& n) { v.resize(n); }
};

}} // namespace cctbx::vector

#endif // CCTBX_VECTOR_MIXED_TRAITS_H
