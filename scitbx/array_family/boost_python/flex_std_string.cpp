#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>

#include <boost/python/overloads.hpp>
#include <boost/python/str.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  // Based on Python-2.4.2/Objects/stringobject.c string_splitlines()
  af::shared<std::string>
  split_lines(
    boost::python::str const& multi_line_string,
    bool keep_ends=false,
    bool count_lines_first=true)
  {
    af::shared<std::string> result;
    PyObject* str_ptr = multi_line_string.ptr();
    const char* data = PyString_AS_STRING(str_ptr);
    std::size_t len = PyString_GET_SIZE(str_ptr);
    std::size_t n_lines = 0;
    for(unsigned i_pass=(count_lines_first ? 0 : 1);i_pass<2;i_pass++) {
      int i = 0;
      int j = 0;
      while (i < len) {
        while (i < len && data[i] != '\n' && data[i] != '\r') i++;
        int eol = i;
        if (i < len) {
          if (data[i] == '\r' && i+1 < len && data[i+1] == '\n') {
            i += 2;
          }
          else {
            i++;
          }
          if (keep_ends) {
            eol = i;
          }
        }
        if (i_pass == 0) {
          n_lines++;
        }
        else {
          result.push_back(std::string(data+j, data+eol));
        }
        j = i;
      }
      if (j < len) {
        if (i_pass == 0) {
          n_lines++;
        }
        else {
          result.push_back(std::string(data+j, data+len));
        }
      }
      if (i_pass == 0) {
        result.reserve(n_lines);
      }
    }
    return result;
  }

  BOOST_PYTHON_FUNCTION_OVERLOADS(split_lines_overloads, split_lines, 1, 3)

} // namespace <anonymous>

  void wrap_flex_std_string()
  {
    using namespace boost::python;
    typedef flex_wrapper<std::string> fw;
    fw::ordered("std_string", scope())
      .def_pickle(flex_pickle_double_buffered<std::string>())
      .def("count", fw::count)
    ;
    def("split_lines", split_lines, split_lines_overloads((
      arg_("multi_line_string"),
      arg_("keep_ends")=false,
      arg_("count_lines_first")=true)));
  }

}}} // namespace scitbx::af::boost_python
