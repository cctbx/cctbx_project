// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_VECTOR_SHARED_STORAGE_TRAITS_H
#define CCTBX_VECTOR_SHARED_STORAGE_TRAITS_H

#include <cctbx/shared_storage.h>

namespace cctbx { namespace vector {

template <class T, class U>
struct algebra_traits<shared_storage<T>, shared_storage<U> >
{
  typedef typename algebra_traits<T, U>::promotion_type value_type;
  typedef shared_storage<value_type> promotion_type;
  typedef shared_storage<bool> bool_type;
  typedef qualifier_v_v type_qualifier;
};

template <class T, class U>
struct algebra_traits<shared_storage<T>, U>
{
  typedef typename algebra_traits<T, U>::promotion_type value_type;
  typedef shared_storage<value_type> promotion_type;
  typedef shared_storage<bool> bool_type;
  typedef qualifier_v_s type_qualifier;
};

template <class T, class U>
struct algebra_traits<T, shared_storage<U> >
{
  typedef typename algebra_traits<T, U>::promotion_type value_type;
  typedef shared_storage<value_type> promotion_type;
  typedef shared_storage<bool> bool_type;
  typedef qualifier_s_v type_qualifier;
};

}} // namespace cctbx::vector

#endif // CCTBX_VECTOR_SHARED_STORAGE_TRAITS_H
