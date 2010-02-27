#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <scitbx/misc/split_lines.h>

#include <boost/python/str.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>

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

  af::shared<std::string>
  strip(
    af::const_ref<std::string> const& self)
  {
    af::shared<std::string> result((reserve(self.size())));
    for (std::size_t i = 0; i < self.size(); i++) {
      std::string trimmed = boost::algorithm::trim_copy(self[i]);
      result.push_back(trimmed);
    }
    SCITBX_ASSERT(result.size() == self.size());
    return result;
  }

  af::shared<std::string>
  upper(
    af::const_ref<std::string> const& self)
  {
    af::shared<std::string> result((reserve(self.size())));
    for (std::size_t i = 0; i <self.size(); i++) {
      std::string uppercase_string = boost::to_upper_copy(self[i]);
      result.push_back(uppercase_string);
    }
    SCITBX_ASSERT(result.size() == self.size());
    return result;
  }

  af::shared<std::string>
  lower(
    af::const_ref<std::string> const& self)
  {
    af::shared<std::string> result((reserve(self.size())));
    for (std::size_t i = 0; i <self.size(); i++) {
      std::string lowercase_string = boost::to_lower_copy(self[i]);
      result.push_back(lowercase_string);
    }
    SCITBX_ASSERT(result.size() == self.size());
    return result;
  }

} // namespace <anonymous>

  void wrap_flex_std_string()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef flex_wrapper<std::string> fw;
    fw::ordered("std_string", scope())
      .def_pickle(flex_pickle_double_buffered<std::string>())
      .def("count", fw::count)
      .def("strip", strip)
      .def("upper", upper)
      .def("lower", lower)
    ;
    def("split_lines", split_lines_wrapper, (
      arg("multi_line_string"),
      arg("keep_ends")=false,
      arg("count_lines_first")=true));
  }

}}} // namespace scitbx::af::boost_python
