/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/ref_c_grid_flex_conversions.h>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_long()
  {
    flex_wrapper<long>::integer("long", boost::python::scope())
      .def_pickle(flex_pickle_single_buffered<long, 21>());

    default_ref_c_grid_flex_conversions<long>();
  }

}}} // namespace scitbx::af::boost_python
