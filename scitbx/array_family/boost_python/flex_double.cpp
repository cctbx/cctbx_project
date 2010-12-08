#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/byte_str.h>
#include <scitbx/array_family/boost_python/range_wrappers.h>
#include <scitbx/array_family/boost_python/numpy_bridge.hpp>
#include <scitbx/math/utils.h>
#include <scitbx/matrix/norms.h>
#include <boost/python/args.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <vector>
#include <set>
#include "flex_helpers.h"

namespace scitbx { namespace af {

namespace {

  /* This will also allow construction from cif-style numeric strings
     where the esd associated with the value is also given in brackets.
     For now the value in the brackets (if present) is ignored
     (although this is checked to be a valid integer).  E.g.
       the string "1.3(2)" will be interpreted as the double 1.3
   */

  flex<double>::type*
  from_std_string(const_ref<std::string> const& strings)
  {
    shared<double> result(reserve(strings.size()));
    for(std::size_t i=0;i<strings.size();i++) {
      std::string s = strings[i];
      std::size_t open_bracket_i = s.find_first_of('(');
      std::size_t close_bracket_i = s.find_last_of(')');
      double value;
      if (open_bracket_i == std::string::npos) {
        value = boost::lexical_cast<double>(s);
      }
      else {
        if (close_bracket_i != s.size()-1) {
          throw std::runtime_error("Unexpected trailing characters after ')'");
        }
        value = boost::lexical_cast<double>(s.substr(0, open_bracket_i));
        // check that value between brackets is a valid integer
        boost::lexical_cast<int>(
          s.substr(open_bracket_i+1, close_bracket_i-open_bracket_i-1));
      }
      result.push_back(value);
    }
    return new flex<double>::type(result, result.size());
  }

  flex<double>::type*
  from_stl_vector_double(std::vector<double> const& v)
  {
    shared<double> result(reserve(v.size()));
    for(std::size_t i=0;i<v.size();i++) {
      result.push_back(v[i]);
    }
    return new flex<double>::type(result, result.size());
  }

  flex<double>::type*
  from_list_or_tuple_of_lists_or_tuples(PyObject* matrix_ptr)
  {
    static const char* argument_error
      = "argument must be a Python list or tuple of lists or tuples.";
    static const char* column_size_error
      = "matrix columns must have identical sizes.";
    shared<double> result;
    std::size_t n_rows = 0;
    std::size_t n_columns = 0;
    if (PyList_Check(matrix_ptr)) {
      n_rows = PyList_GET_SIZE(matrix_ptr);
      for(std::size_t i_row=0;i_row<n_rows;i_row++) {
        PyObject* row = PyList_GET_ITEM(matrix_ptr, i_row);
        if (PyList_Check(row)) {
          if (i_row == 0) {
            n_columns = PyList_GET_SIZE(row);
            result.reserve(n_rows*n_columns);
          }
          else {
            if (PyList_GET_SIZE(row) != n_columns) {
              throw std::runtime_error(column_size_error);
            }
          }
          for(std::size_t i_column=0;i_column<n_columns;i_column++) {
            PyObject* elem = PyList_GET_ITEM(row, i_column);
            double value = PyFloat_AsDouble(elem);
            if (PyErr_Occurred()) boost::python::throw_error_already_set();
            result.push_back(value);
          }
        }
        else if (PyTuple_Check(row)) {
          if (i_row == 0) {
            n_columns = PyTuple_GET_SIZE(row);
          }
          else {
            if (PyTuple_GET_SIZE(row) != n_columns) {
              throw std::runtime_error(column_size_error);
            }
          }
          for(std::size_t i_column=0;i_column<n_columns;i_column++) {
            PyObject* elem = PyTuple_GET_ITEM(row, i_column);
            double value = PyFloat_AsDouble(elem);
            if (PyErr_Occurred()) boost::python::throw_error_already_set();
            result.push_back(value);
          }
        }
        else if (i_row == 0) {
          result.reserve(n_rows);
          for (;i_row<n_rows;i_row++) {
            PyObject* elem = PyList_GET_ITEM(matrix_ptr, i_row);
            double value = PyFloat_AsDouble(elem);
            if (PyErr_Occurred()) boost::python::throw_error_already_set();
            result.push_back(value);
          }
          return new flex<double>::type(result, flex_grid<>(n_rows));
        }
        else {
          throw std::runtime_error(argument_error);
        }
      }
    }
    else if (PyTuple_Check(matrix_ptr)) {
      n_rows = PyTuple_GET_SIZE(matrix_ptr);
      for(std::size_t i_row=0;i_row<n_rows;i_row++) {
        PyObject* row = PyTuple_GET_ITEM(matrix_ptr, i_row);
        if (PyList_Check(row)) {
          if (i_row == 0) {
            n_columns = PyList_GET_SIZE(row);
            result.reserve(n_rows*n_columns);
          }
          else {
            if (PyList_GET_SIZE(row) != n_columns) {
              throw std::runtime_error(column_size_error);
            }
          }
          for(std::size_t i_column=0;i_column<n_columns;i_column++) {
            PyObject* elem = PyList_GET_ITEM(row, i_column);
            double value = PyFloat_AsDouble(elem);
            if (PyErr_Occurred()) boost::python::throw_error_already_set();
            result.push_back(value);
          }
        }
        else if (PyTuple_Check(row)) {
          if (i_row == 0) {
            n_columns = PyTuple_GET_SIZE(row);
          }
          else {
            if (PyTuple_GET_SIZE(row) != n_columns) {
              throw std::runtime_error(column_size_error);
            }
          }
          for(std::size_t i_column=0;i_column<n_columns;i_column++) {
            PyObject* elem = PyTuple_GET_ITEM(row, i_column);
            double value = PyFloat_AsDouble(elem);
            if (PyErr_Occurred()) boost::python::throw_error_already_set();
            result.push_back(value);
          }
        }
        else if (i_row == 0) {
          result.reserve(n_rows);
          for (;i_row<n_rows;i_row++) {
            PyObject* elem = PyTuple_GET_ITEM(matrix_ptr, i_row);
            double value = PyFloat_AsDouble(elem);
            if (PyErr_Occurred()) boost::python::throw_error_already_set();
            result.push_back(value);
          }
          return new flex<double>::type(result, flex_grid<>(n_rows));
        }
        else {
          throw std::runtime_error(argument_error);
        }
      }
    }
    else {
      throw std::runtime_error(argument_error);
    }
    if (n_rows == 0) {
      return new flex<double>::type(result, flex_grid<>(0));
    }
    return new flex<double>::type(result, flex_grid<>(n_rows, n_columns));
  }

  flex<double>::type*
  from_list_of_lists_or_tuples(boost::python::list const& matrix)
  {
    return from_list_or_tuple_of_lists_or_tuples(matrix.ptr());
  }

  flex<double>::type*
  from_tuple_of_lists_or_tuples(boost::python::tuple const& matrix)
  {
    return from_list_or_tuple_of_lists_or_tuples(matrix.ptr());
  }

  shared<double>
  extract_double_attributes(
    boost::python::object array,
    const char* attribute_name,
    boost::python::object none_substitute)
  {
    PyObject* array_ptr = array.ptr();
#if PY_VERSION_HEX >= 0x02050000
    const char* attr_name = attribute_name;
#else
    char* attr_name = const_cast<char*>(attribute_name);
#endif
    PyObject* none_substitute_ptr = none_substitute.ptr();
    if (PyList_Check(array_ptr)) {
      std::size_t len_array = PyList_GET_SIZE(array_ptr);
      shared<double> result((reserve(len_array)));
      for(std::size_t i=0;i<len_array;i++) {
        PyObject* elem = PyList_GET_ITEM(array_ptr, i);
        PyObject* elem_attr = PyObject_GetAttrString(elem, attr_name);
        if (!elem_attr) boost::python::throw_error_already_set();
        if (elem_attr == Py_None) elem_attr = none_substitute_ptr;
        double value = PyFloat_AsDouble(elem_attr);
        if (PyErr_Occurred()) boost::python::throw_error_already_set();
        result.push_back(value);
      }
      return result;
    }
    if (PyTuple_Check(array_ptr)) {
      std::size_t len_array = PyTuple_GET_SIZE(array_ptr);
      shared<double> result((reserve(len_array)));
      for(std::size_t i=0;i<len_array;i++) {
        PyObject* elem = PyTuple_GET_ITEM(array_ptr, i);
        PyObject* elem_attr = PyObject_GetAttrString(elem, attr_name);
        if (!elem_attr) boost::python::throw_error_already_set();
        if (elem_attr == Py_None) elem_attr = none_substitute_ptr;
        double value = PyFloat_AsDouble(elem_attr);
        if (PyErr_Occurred()) boost::python::throw_error_already_set();
        result.push_back(value);
      }
      return result;
    }
    throw error("array must be a Python list or tuple.");
  }

  versa<std::complex<double>, flex_grid<> >
  mul_ar_sc(
    versa<double, flex_grid<> > const& self,
    std::complex<double> const& rhs)
  {
    versa<std::complex<double>, flex_grid<> > result(
      self.accessor(),
      init_functor_null<std::complex<double> >());
    std::complex<double>* r = result.begin();
    const double* s = self.begin();
    const double* s_end = self.end();
    while (s != s_end) {
      *r++ = (*s++) * rhs;
    }
    return result;
  }

  bool
  all_approx_equal_a_a(
    const_ref<double> const& self,
    const_ref<double> const& other,
    double tolerance=1e-6)
  {
    return self.all_approx_equal(other, tolerance);
  }

  bool
  all_approx_equal_a_s(
    const_ref<double> const& self,
    double other,
    double tolerance=1e-6)
  {
    return self.all_approx_equal(other, tolerance);
  }

  bool
  all_approx_equal_relatively_a_a(
    const_ref<double> const& self,
    const_ref<double> const& other,
    double relative_error=1e-6)
  {
    return self.all_approx_equal_relatively(other, relative_error);
  }

  bool
  all_approx_equal_relatively_a_s(
    const_ref<double> const& self,
    double other,
    double relative_error=1e-6)
  {
    return self.all_approx_equal_relatively(other, relative_error);
  }

  af::versa<float, af::flex_grid<> >
  as_float(
    af::const_ref<double, af::flex_grid<> > const& O)
  {
    af::versa<float, af::flex_grid<> > result(
      O.accessor(), af::init_functor_null<float>());
    std::size_t n = O.accessor().size_1d();
    float* r = result.begin();
    for(std::size_t i=0;i<n;i++) {
      r[i] = static_cast<float>(O[i]);
    }
    return result;
  }

  double norm_1_a(const_ref<double> const &self) {
    return matrix::norm_1(self);
  }

  /* For allowed syntax for the optional format_string argument see:
       http://www.boost.org/libs/format/doc/format.html#syntax
   */
  af::shared<std::string>
  as_string(af::const_ref<double, af::flex_grid<> > const& O,
            std::string format_string="%d")
  {
    af::shared<std::string> result((reserve(O.size())));
    std::size_t n = O.accessor().size_1d();
    for(std::size_t i=0;i<n;i++) {
      result.push_back((boost::format(format_string) %O[i]).str());
    }
    return result;
  }

  shared<double>
  round(
    const_ref<double> const& self,
    int n_digits=0)
  {
    shared<double> result(self.size(), init_functor_null<double>());
    for(std::size_t i=0;i<self.size();i++) {
      result[i] = math::round(self[i], n_digits);
    }
    return result;
  }

  template <typename S>
  shared<double>
  select_stl_iterable(
    versa<double, flex_grid<> > const& self,
    S const& selection)
  {
    shared<double> result(selection.size(), init_functor_null<double>());
    typename S::const_iterator sel_end = selection.end();
    double* r = result.begin();
    const double* s = self.begin();
    typedef typename S::value_type svt;
    svt self_size = boost::numeric_cast<svt>(self.size());
    for(typename S::const_iterator sel=selection.begin();sel!=sel_end;sel++) {
      SCITBX_ASSERT(*sel < self_size);
      *r++ = s[*sel];
    }
    SCITBX_ASSERT(r == result.end());
    return result;
  }

  std::string mathematica_form(af::const_ref<double, flex_grid<> > const &self)
  {
    /* Let's use Mathematica cleverness instead of working it out in C++ */
    std::ostringstream o, c;
    c << "{";
    for (std::size_t i=0; i<self.size(); ++i) {
      c << self[i];
      if (i != self.size() - 1) c << ",";
    }
    c << "}";
    std::string coeffs = c.str();
    boost::replace_all(coeffs, "e", "*^");
    if (self.accessor().nd() > 1) o << "Fold[Partition,";
    o << coeffs;
    if (self.accessor().nd() > 1) {
      o << ",";
      if (self.accessor().nd() > 2) o << "Reverse[";
      o << "{";
      flex_grid<>::index_type indices = self.accessor().all();
      for (int i=1; i<indices.size(); ++i) {
        o << indices[i];
        if (i != indices.size() - 1) o << ",";
      }
      o << "}";
      if (self.accessor().nd() > 2) o << "]";
      o << "]";
    }
    return o.str();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_flex_double_matrix(
    flex_wrapper<double>::class_f_t& class_f_t);

  void wrap_flex_double()
  {
    using namespace boost::python;
    using boost::python::arg;

    typedef flex_wrapper<double> f_w;
    f_w::class_f_t class_f_t(f_w::numeric("double", scope()));
    class_f_t
      .def_pickle(flex_pickle_single_buffered<double>())
      .def("copy_to_byte_str", copy_to_byte_str<versa<double, flex_grid<> > >)
      .def("slice_to_byte_str",
        slice_to_byte_str<versa<double,flex_grid<> > >)
      .def("__init__", make_constructor(
        from_std_string, default_call_policies()))
      .def("__init__", make_constructor(
        from_stl_vector_double, default_call_policies()))
      .def("__init__", make_constructor(
        from_list_of_lists_or_tuples, default_call_policies()))
      .def("__init__", make_constructor(
        from_tuple_of_lists_or_tuples, default_call_policies()))
      .def("__init__", make_constructor(
        flex_double_from_numpy_array, default_call_policies()))
      .def("__mul__", mul_ar_sc)
      .def("__rmul__", mul_ar_sc)
      .def("__mul__", f_w::mul_a_s) // re-define so it is found first
      .def("__rmul__", f_w::mul_a_s) // re-define so it is found first
      .def("add_selected", add_selected_bool_a<double>, (
        arg("flags"), arg("values")))
      .def("add_selected", add_selected_bool_s<double>, (
        arg("flags"), arg("value")))
      .def("add_selected", add_selected_unsigned_a<double, std::size_t>, (
        arg("indices"), arg("values")))
      .def("add_selected", add_selected_unsigned_s<double, std::size_t>, (
        arg("indices"), arg("value")))
      .def("all_approx_equal",
        all_approx_equal_a_a, (
          arg("other"),
          arg("tolerance")=1e-6))
      .def("all_approx_equal",
        all_approx_equal_a_s, (
          arg("other"),
          arg("tolerance")=1e-6))
      .def("all_approx_equal_relatively",
        all_approx_equal_relatively_a_a, (
          arg("other"),
          arg("relative_error")=1e-6))
      .def("all_approx_equal_relatively",
        all_approx_equal_relatively_a_s, (
          arg("other"),
          arg("relative_error")=1e-6))
      .def("as_float", as_float)
      .def("as_string", as_string, (
          arg("other"),
          arg("format_string")="%d"))
      .def("mathematica_form", mathematica_form)
      .def("round", round, (arg("n_digits")=0))
      .def("select", select_stl_iterable<std::vector<unsigned> >, (
        arg("selection")))
      .def("select", select_stl_iterable<std::set<unsigned> >, (
        arg("selection")))
      .def("norm_1", norm_1_a)
      .def("as_numpy_array", flex_double_as_numpy_array, (
        arg("optional")=false))
    ;
    def(
      "double_from_byte_str",
      shared_from_byte_str<double>,
      (arg("byte_str")));
    range_wrappers<double, long>::wrap("double_range");

    typedef return_value_policy<return_by_value> rbv;
    typedef af::min_max_mean<double> mmm;
    class_<mmm>("min_max_mean_double", no_init)
      .def(init<af::const_ref<double> const&>((arg("values"))))
      .def_readonly("n", &mmm::n)
      .add_property("min", make_getter(&mmm::min, rbv()))
      .add_property("max", make_getter(&mmm::max, rbv()))
      .add_property("sum", make_getter(&mmm::sum, rbv()))
      .add_property("mean", make_getter(&mmm::mean, rbv()))
    ;

    def("extract_double_attributes", extract_double_attributes,
      (arg("array"), arg("attribute_name"), arg("none_substitute")));

    wrap_flex_double_matrix(class_f_t);
  }

}}} // namespace scitbx::af::boost_python
