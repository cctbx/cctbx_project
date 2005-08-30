#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <cctbx/sgtbx/lattice_symmetry.h>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_evaluated_axis()
  {
    flex_wrapper<cctbx::sgtbx::lattice_symmetry::evaluated_axis_t>::plain(
     "evaluated_axis_t");
  }

}}} // namespace scitbx::af::boost_python
