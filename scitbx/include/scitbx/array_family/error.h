/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Mar: moved parts from misc.h (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_ERROR_H
#define SCITBX_ARRAY_FAMILY_ERROR_H

#include <stdexcept>

// FIXES for broken compilers
#include <boost/config.hpp>

namespace scitbx { namespace af {

  inline
  void throw_range_error() {
    throw std::range_error("array_family");
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ERROR_H
