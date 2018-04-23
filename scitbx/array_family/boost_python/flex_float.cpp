#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/range_wrappers.h>
#include <boost/python/make_constructor.hpp>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_float()
  {
    using namespace boost::python;
    flex_wrapper<float>::numeric("float", boost::python::scope())
      .def_pickle(flex_pickle_single_buffered<float>())
    ;
    range_wrappers<float, int>::wrap("float_range");
  }

}}} // namespace scitbx::af::boost_python
