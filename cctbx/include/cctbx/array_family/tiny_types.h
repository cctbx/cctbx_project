// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_TINY_TYPES_H
#define CCTBX_ARRAY_FAMILY_TINY_TYPES_H

#include <cctbx/array_family/tiny.h>

namespace cctbx { namespace af {

  typedef tiny<int, 3> int3;
  typedef tiny<int, 9> int9;
  typedef tiny<double, 2> double2;
  typedef tiny<double, 3> double3;
  typedef tiny<double, 6> double6;
  typedef tiny<double, 9> double9;

}} // namespace cctbx::af

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace boost {

  template<>
  struct has_trivial_destructor<cctbx::af::int3> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<cctbx::af::int9> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<cctbx::af::double2> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<cctbx::af::double3> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<cctbx::af::double6> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<cctbx::af::double9> {
    static const bool value = true;
  };

}

#endif // !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#endif // CCTBX_ARRAY_FAMILY_TINY_TYPES_H
