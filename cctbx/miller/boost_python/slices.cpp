
#include <cctbx/miller/slices.h>
#include <boost/python/def.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct slice_wrappers
  {

    static void
    wrap()
    {
      using namespace boost::python;
      def("simple_slice", simple_slice, (
        arg("indices"),
        arg("slice_axis"),
        arg("slice_index")));
      def("multi_slice", multi_slice, (
        arg("indices"),
        arg("slice_axis"),
        arg("slice_start"),
        arg("slice_end")));
    }
  };

} // namespace <anoymous>

  void wrap_slices()
  {
    slice_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
