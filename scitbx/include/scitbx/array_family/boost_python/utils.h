/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/bpl_utils.h (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_UTILS_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_UTILS_H

#include <scitbx/array_family/flex_grid_accessor.h>
#include <scitbx/array_family/versa.h>

namespace scitbx { namespace af { namespace boost_python {

  void raise_must_be_0_based_1d();
  void raise_must_be_0_based_3d();
  void assert_0_based_1d(flex_grid<> const& grid);
  void assert_0_based_3d(flex_grid<> const& grid);
  void raise_shared_size_mismatch();
  void raise_incompatible_arrays();

  template <typename ElementType>
  shared_plain<ElementType>
  flex_as_base_array(versa<ElementType, flex_grid<> > const& a)
  {
    if (!a.check_shared_size()) raise_shared_size_mismatch();
    assert_0_based_1d(a.accessor());
    shared_plain<ElementType> b = a.as_base_array();
    if (a.size() != b.size()) raise_shared_size_mismatch();
    return b;
  }

}}} // namespace scitbx::af::boost_python

#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_UTILS_H
