// $Id$
/* Copyright (c) 2002 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_TYPE_TRAITS_H
#define CCTBX_ARRAY_FAMILY_TYPE_TRAITS_H

#include <boost/config.hpp>

#if defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#define CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR

#else

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#define CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR \
{ \
  BOOST_STATIC_ASSERT(::boost::has_trivial_destructor<ElementType>::value); \
}

#include <complex>

namespace boost {
  template <typename T>
  struct has_trivial_destructor<std::complex<T> > {
    // we really hope that this is true ...
    static const bool value = true;
  };
}

#endif // defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#endif // CCTBX_ARRAY_FAMILY_TYPE_TRAITS_H
