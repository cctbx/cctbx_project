// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: moved parts from misc.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_ERROR_H
#define CCTBX_ARRAY_FAMILY_ERROR_H

#include <stdexcept>

// FIXES for broken compilers
#include <boost/config.hpp>

namespace cctbx { namespace af {

  inline
  void throw_range_error() {
    throw std::range_error("array_family");
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_ERROR_H
