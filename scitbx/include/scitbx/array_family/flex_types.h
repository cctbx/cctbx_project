/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_FLEX_TYPES_H
#define SCITBX_ARRAY_FAMILY_FLEX_TYPES_H

#include <complex>
#include <scitbx/array_family/flex_grid_accessor.h>
#include <scitbx/array_family/versa.h>

namespace scitbx { namespace af {

  typedef versa<bool, flex_grid<> > flex_bool;
  typedef versa<int, flex_grid<> > flex_int;
  typedef versa<long, flex_grid<> > flex_long;
  typedef versa<std::size_t, flex_grid<> > flex_size_t;
  typedef versa<float, flex_grid<> > flex_float;
  typedef versa<double, flex_grid<> > flex_double;
  typedef versa<std::complex<double>, flex_grid<> > flex_complex_double;

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_FLEX_TYPES_H
