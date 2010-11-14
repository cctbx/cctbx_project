#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/range_wrappers.h>
#include <scitbx/array_family/boost_python/numpy_bridge.hpp>
#include <boost/python/make_constructor.hpp>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_float()
  {
    using namespace boost::python;
    flex_wrapper<float>::numeric("float", boost::python::scope())
      .def_pickle(flex_pickle_single_buffered<float>())
      .def("__init__", make_constructor(
        flex_float_from_numpy_array, default_call_policies()))
      .def("as_numpy_array", flex_float_as_numpy_array, (
        boost::python::arg("optional")=false))
    ;
    range_wrappers<float, int>::wrap("float_range");
  }

}}} // namespace scitbx::af::boost_python
