#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <scitbx/misc/split_lines.h>

#include <boost/python/dict.hpp>
#include <boost/python/str.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>

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

  boost::python::dict
  i_seqs_by_value(
    af::const_ref<std::string> const& self)
  {
    boost::python::dict result;
    typedef boost::unordered_map<std::string, std::size_t> map_s_s_t;
    map_s_s_t counts;
    for(std::size_t i_seq=0;i_seq<self.size();i_seq++) {
      counts[self[i_seq]]++;
    }
    boost::scoped_array<af::shared<std::size_t> > buffer(
      new af::shared<std::size_t>[counts.size()]);
    map_s_s_t buffer_indices;
    {
      std::size_t ib = 0;
      typedef map_s_s_t::const_iterator it;
      it e = counts.end();
      for(it i=counts.begin();i!=e;i++,ib++) {
        buffer[ib].reserve(i->second);
        buffer_indices[i->first] = ib;
      }
      SCITBX_ASSERT(buffer_indices.size() == counts.size());
    }
    for(std::size_t i_seq=0;i_seq<self.size();i_seq++) {
      std::size_t ib = buffer_indices[self[i_seq]];
      buffer[ib].push_back(i_seq);
    }
    {
      typedef map_s_s_t::const_iterator it;
      it e = buffer_indices.end();
      for(it i=buffer_indices.begin();i!=e;i++) {
        result[i->first] = buffer[i->second];
      }
    }
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
      .def("i_seqs_by_value", i_seqs_by_value)
    ;
    def("split_lines", split_lines_wrapper, (
      arg("multi_line_string"),
      arg("keep_ends")=false,
      arg("count_lines_first")=true));
  }

}}} // namespace scitbx::af::boost_python
