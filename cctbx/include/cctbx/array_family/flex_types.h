// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_FLEX_TYPES_H
#define CCTBX_ARRAY_FAMILY_FLEX_TYPES_H

#include <cctbx/array_family/flex_grid_accessor.h>
#include <cctbx/array_family/versa.h>

namespace cctbx { namespace af {

  typedef versa<bool, flex_grid<> > flex_bool;
  typedef versa<int, flex_grid<> > flex_int;
  typedef versa<long, flex_grid<> > flex_long;
  typedef versa<std::size_t, flex_grid<> > flex_size_t;
  typedef versa<float, flex_grid<> > flex_float;
  typedef versa<double, flex_grid<> > flex_double;
  typedef versa<std::complex<double>, flex_grid<> > flex_complex_double;

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_FLEX_TYPES_H
