#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/versa_matrix.h>
#include <boost/python/args.hpp>
#include <boost/python/make_constructor.hpp>

namespace scitbx { namespace af {

namespace {

  flex<double>::type*
  from_stl_vector_double(std::vector<double> const& v)
  {
    af::shared<double> result(af::reserve(v.size()));
    for(std::size_t i=0;i<v.size();i++) {
      result.push_back(v[i]);
    }
    return new flex<double>::type(result, result.size());
  }

  af::shared<double>
  extract_double_attributes(
    boost::python::object array,
    const char* attribute_name,
    boost::python::object none_substitute)
  {
    PyObject* array_ptr = array.ptr();
    char* attr_name = const_cast<char*>(attribute_name);
    PyObject* none_substitute_ptr = none_substitute.ptr();
    if (PyList_Check(array_ptr)) {
      std::size_t len_array = PyList_GET_SIZE(array_ptr);
      af::shared<double> result((reserve(len_array)));
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
      af::shared<double> result((reserve(len_array)));
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

} // namespace <anonymous>

namespace boost_python {

  void wrap_flex_double()
  {
    using namespace boost::python;

    flex_wrapper<double>::numeric("double", scope())
      .def_pickle(flex_pickle_single_buffered<double>())
      .def("__init__", make_constructor(
        from_stl_vector_double, default_call_policies()))
      .def("matrix_diagonal",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&)) matrix_diagonal)
      .def("matrix_diagonal_sum",
        (double(*)(
          const_ref<double, c_grid<2> > const&)) matrix_diagonal_sum)
      .def("matrix_trace",
        (double(*)(
          const_ref<double, c_grid<2> > const&)) matrix_diagonal_sum)
      .def("matrix_diagonal_product",
        (double(*)(
          const_ref<double, c_grid<2> > const&)) matrix_diagonal_product)
      .def("matrix_multiply",
        (versa<double, c_grid<2> >(*)(
          const_ref<double, c_grid<2> > const&,
          const_ref<double, c_grid<2> > const&)) matrix_multiply)
      .def("matrix_multiply",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&,
          const_ref<double> const&)) matrix_multiply)
      .def("matrix_multiply",
        (shared<double>(*)(
          const_ref<double> const&,
          const_ref<double, c_grid<2> > const&)) matrix_multiply)
      .def("matrix_multiply",
        (double(*)(
          const_ref<double> const&,
          const_ref<double> const&)) matrix_multiply)
      .def("dot",
        (double(*)(
          const_ref<double> const&,
          const_ref<double> const&)) matrix_multiply)
      .def("transpose_in_place",
        (void(*)(versa<double, flex_grid<> >&)) transpose_in_place)
      .def("matrix_lu_decomposition_in_place",
        (shared<std::size_t>(*)(
          ref<double, c_grid<2> > const&)) matrix_lu_decomposition_in_place)
      .def("matrix_lu_back_substitution",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&,
          const_ref<std::size_t> const&,
          const_ref<double> const&)) matrix_lu_back_substitution, (
        arg_("pivot_indices"), arg_("b")))
    ;

    def("extract_double_attributes", extract_double_attributes,
      (arg_("array"), arg_("attribute_name"), arg_("none_substitute")));
  }

}}} // namespace scitbx::af::boost_python
