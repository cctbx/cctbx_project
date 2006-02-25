#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <scitbx/misc/split_lines.h>

#include <boost/python/overloads.hpp>
#include <boost/python/str.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  af::shared<std::string>
  split_lines_wrapper(
    boost::python::str const& multi_line_string,
    bool keep_ends=false,
    bool count_lines_first=true)
  {
    PyObject* str_ptr = multi_line_string.ptr();
    return misc::split_lines(
      PyString_AS_STRING(str_ptr),
      PyString_GET_SIZE(str_ptr),
      keep_ends,
      count_lines_first);
  }

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    split_lines_overloads, split_lines_wrapper, 1, 3)

} // namespace <anonymous>

  void wrap_flex_std_string()
  {
    using namespace boost::python;
    typedef flex_wrapper<std::string> fw;
    fw::ordered("std_string", scope())
      .def_pickle(flex_pickle_double_buffered<std::string>())
      .def("count", fw::count)
    ;
    def("split_lines", split_lines_wrapper, split_lines_overloads((
      arg_("multi_line_string"),
      arg_("keep_ends")=false,
      arg_("count_lines_first")=true)));
  }

}}} // namespace scitbx::af::boost_python
