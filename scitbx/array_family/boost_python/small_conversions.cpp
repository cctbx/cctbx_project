/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Fragments from cctbx/misc/bpl_utils.cpp (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/boost_python/small_conversions.h>
#include <scitbx/array_family/flex_grid_accessor.h>

namespace scitbx { namespace af { namespace boost_python {

  void register_small_conversions()
  {
    static bool already_registered = false;
    if (already_registered) return;
    already_registered = true;

    scitbx::boost_python::container_conversions::tuple_mapping_fixed_capacity<
      flex_grid_default_index_type>();
  }

}}} // namespace scitbx::af::boost_python
