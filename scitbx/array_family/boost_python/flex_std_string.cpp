#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_std_string()
  {
    typedef flex_wrapper<std::string> fw;
    fw::ordered("std_string", boost::python::scope())
      .def_pickle(flex_pickle_double_buffered<std::string>())
      .def("count", fw::count)
    ;
  }

}}} // namespace scitbx::af::boost_python
