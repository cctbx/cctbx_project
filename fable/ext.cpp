#include <boost/python.hpp>

#include <fem/format.hpp>
#include <fem/utils/equivalence.hpp>
#include <fem/utils/string_to_double_fmt.hpp>
#include <boost_adaptbx/error_utils.h>

#include <cctype>

namespace fable { namespace ext {

  namespace bp = boost::python;

  // helper, not available in Python
  // compare with fable/__init__.py
  unsigned
  get_code_stop(
    bp::object const& code,
    int stop)
  {
    unsigned len_code = bp::len(code);
    if (stop < 0) return len_code;
    ASSERTBX(stop <= len_code);
    return static_cast<unsigned>(stop);
  }

  // compare with fable/__init__.py
  int
  unsigned_integer_scan(
    bp::object const& code,
    unsigned start,
    int stop)
  {
    unsigned code_stop = get_code_stop(code, stop);
    char const* s = bp::extract<char const*>(code)();
    return fem::utils::unsigned_integer_scan(s, start, code_stop);
  }

  // compare with fable/__init__.py
  int
  floating_point_scan_after_exponent_char(
    bp::object const& code,
    unsigned start,
    int stop)
  {
    unsigned code_stop = get_code_stop(code, stop);
    char const* s = bp::extract<char const*>(code)();
    unsigned i = start;
    if (i < code_stop) {
      int c = s[i];
      if (c == '+' || c == '-') {
        i += 1;
      }
      return unsigned_integer_scan(code, i, stop);
    }
    return -1;
  }

  // compare with fable/__init__.py
  int
  floating_point_scan_after_dot(
    bp::object const& code,
    unsigned start,
    int stop)
  {
    unsigned code_stop = get_code_stop(code, stop);
    char const* s = bp::extract<char const*>(code)();
    int i = unsigned_integer_scan(code, start, stop);
    if (i < 0) i = start;
    if (i < code_stop) {
      int c = s[i];
      if (c == 'e' || c == 'd') {
        return floating_point_scan_after_exponent_char(code, i+1, stop);
      }
    }
    return i;
  }

  // compare with fable/__init__.py
  int
  identifier_scan(
    bp::object const& code,
    unsigned start,
    int stop)
  {
    unsigned code_stop = get_code_stop(code, stop);
    char const* s = bp::extract<char const*>(code)();
    unsigned i = start;
    if (i < code_stop) {
      int c = s[i++];
      if ((c < 'a' || c > 'z') && c != '_') return -1;
      while (i < code_stop) {
        c = s[i++];
        if (   (c < 'a' || c > 'z')
            && (c < '0' || c > '9') && c != '_') return i-1;
      }
      return i;
    }
    return -1;
  }

  // compare with fable/__init__.py
  int
  find_closing_parenthesis(
    bp::object const& code,
    unsigned start,
    int stop)
  {
    unsigned code_stop = get_code_stop(code, stop);
    char const* s = bp::extract<char const*>(code)();
    unsigned n_inner = 0;
    for(unsigned i=start;i<code_stop;i++) {
      int c = s[i];
      if (c == ')') {
        if (n_inner == 0) return i;
        n_inner--;
      }
      else if (c == '(') {
        n_inner++;
      }
    }
    return -1;
  }

  bp::list
  exercise_fem_format_tokenizer(
    std::string const& fmt)
  {
    bp::list result;
    fem::format::tokenizer tz(fmt.c_str(), fmt.size());
    typedef std::vector<fem::utils::token>::const_iterator it;
    it e = tz.tokens.end();
    for(it i=tz.tokens.begin();i!=e;i++) {
      result.append(bp::make_tuple(i->type, i->value));
    }
    return result;
  }

  bp::tuple
  exercise_fem_utils_string_to_double(
    std::string const& str)
  {
    fem::utils::simple_istream_from_std_string inp(str.c_str());
    fem::utils::string_to_double conv(inp);
    return bp::make_tuple(
      conv.result,
      (conv.error_message ? bp::object(*conv.error_message) : bp::object()),
      inp.get());
  }

  bp::tuple
  exercise_fem_utils_string_to_double_fmt(
    std::string const& str,
    int w,
    int d,
    bool blanks_zero,
    int exp_scale)
  {
    fem::utils::simple_istream_from_std_string inp(str.c_str());
    fem::utils::string_to_double_fmt conv(inp, w, d, blanks_zero, exp_scale);
    return bp::make_tuple(
      conv.result,
      (conv.error_message ? bp::object(*conv.error_message) : bp::object()),
      inp.get());
  }

  struct equivalence_array_alignment_wrappers
  {
    typedef fem::utils::equivalence::array_alignment w_t;

    static
    bp::list
    infer_diffs0_from_diff_matrix(
      w_t& O)
    {
      O.infer_diffs0_from_diff_matrix();
      bp::list result;
      size_t n = O.diffs0.size();
      for(size_t i=0;i<n;i++) {
        result.append(O.diffs0[i]);
      }
      return result;
    }

    static
    void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("fem_utils_equivalence_array_alignment", no_init)
        .def(init<size_t>((arg("members_size"))))
        .def("add_anchor", &w_t::add_anchor, (
          arg("i0"), arg("a0"), arg("i1"), arg("a1")))
        .def("infer_diffs0_from_diff_matrix", infer_diffs0_from_diff_matrix)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;
    def("unsigned_integer_scan", unsigned_integer_scan, (
      arg("code"), arg("start")=0, arg("stop")=-1));
    def("floating_point_scan_after_exponent_char",
      floating_point_scan_after_exponent_char, (
        arg("code"), arg("start")=0, arg("stop")=-1));
    def("floating_point_scan_after_dot", floating_point_scan_after_dot, (
      arg("code"), arg("start")=0, arg("stop")=-1));
    def("identifier_scan", identifier_scan, (
      arg("code"), arg("start")=0, arg("stop")=-1));
    def("find_closing_parenthesis", find_closing_parenthesis, (
      arg("code"), arg("start")=0, arg("stop")=-1));

    def("exercise_fem_format_tokenizer", exercise_fem_format_tokenizer, (
      arg("fmt")));
    def("exercise_fem_utils_string_to_double",
      exercise_fem_utils_string_to_double, (
        arg("str")));
    def("exercise_fem_utils_string_to_double_fmt",
      exercise_fem_utils_string_to_double_fmt, (
        arg("str"), arg("w"), arg("d"), arg("blanks_zero"), arg("exp_scale")));

    equivalence_array_alignment_wrappers::wrap();
  }

}} // namespace fable::ext

BOOST_PYTHON_MODULE(fable_ext)
{
  fable::ext::init_module();
}
