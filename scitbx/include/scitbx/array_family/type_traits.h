/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Feb: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_TYPE_TRAITS_H
#define SCITBX_ARRAY_FAMILY_TYPE_TRAITS_H

namespace scitbx { namespace af {

  struct false_type {};
  struct true_type {};

}}

#include <boost/config.hpp>

#if defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#define SCITBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR

namespace scitbx { namespace af {

  template <typename T>
  struct has_trivial_destructor {
    typedef false_type value;
  };

}}

#else

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#define SCITBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR \
{ \
  BOOST_STATIC_ASSERT(::boost::has_trivial_destructor<ElementType>::value); \
}

#include <complex>

namespace boost {
  template <typename T>
  struct has_trivial_destructor<std::complex<T> > {
    // we really hope that this is true ...
    static const bool value = ::boost::has_trivial_destructor<T>::value;
  };
}

namespace scitbx { namespace af {

  template <bool T>
  struct bool_as_type;

  template <>
  struct bool_as_type<false> {
    typedef false_type value;
  };

  template <>
  struct bool_as_type<true> {
    typedef true_type value;
  };

  template <typename T>
  struct has_trivial_destructor {
    typedef typename bool_as_type<
      ::boost::has_trivial_destructor<T>::value>::value value;
  };

}}

#endif // defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#endif // SCITBX_ARRAY_FAMILY_TYPE_TRAITS_H
