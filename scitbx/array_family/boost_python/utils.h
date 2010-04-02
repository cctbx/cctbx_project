#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_UTILS_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_UTILS_H

#include <boost/python/detail/wrap_python.hpp>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>

namespace scitbx { namespace af { namespace boost_python {

  void raise_must_be_0_based_1d();
  void raise_must_be_0_based_2d();
  void raise_must_be_0_based_3d();
  void assert_0_based_1d(flex_grid<> const& grid);
  void assert_0_based_2d(flex_grid<> const& grid);
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
