// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: moved parts from ref.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_MISC_H
#define CCTBX_ARRAY_FAMILY_MISC_H

#include <stdexcept>

// FIXES for broken compilers
#include <boost/config.hpp>

// XXX
#include <iostream>
#define CheckPoint std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl << std::flush

namespace cctbx { namespace af {

  inline
  void throw_range_error() {
    throw std::range_error("array_family");
  }

  // XXX use std::copy if compiler permits
  template <typename InputElementType,
            typename OutputElementType>
  OutputElementType*
  copy_typeconv(
    const InputElementType* first,
    const InputElementType* last,
    OutputElementType* result)
  {
    OutputElementType* p = result;
    while (first != last) *p++ = OutputElementType(*first++);
    return result;
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_MISC_H
