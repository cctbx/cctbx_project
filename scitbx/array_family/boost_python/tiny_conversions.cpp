/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Fragments from cctbx/misc/bpl_utils.cpp (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/boost_python/tiny_conversions.h>

namespace scitbx { namespace af { namespace boost_python {

  void register_tiny_types_conversions()
  {
    static bool already_registered = false;
    if (already_registered) return;
    already_registered = true;

    using namespace scitbx::boost_python::container_conversions;
    tuple_mapping_fixed_size<int3>();
    tuple_mapping_fixed_size<int9>();
    tuple_mapping_fixed_size<long3>();
    tuple_mapping_fixed_size<double2>();
    tuple_mapping_fixed_size<double3>();
    tuple_mapping_fixed_size<double6>();
    tuple_mapping_fixed_size<double9>();
  }

}}} // namespace scitbx::af::boost_python
