#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>

#include <scitbx/sparse/matrix.h>

namespace scitbx { namespace sparse { namespace boost_python {

template<typename T>
struct matrix_wrapper
{
  typedef matrix<T> wt;
  typedef typename wt::row_index row_index;
  typedef typename wt::column_index column_index;
  typedef typename wt::value_type value_type;

  static void setitem(wt& self, boost::python::tuple ij, T x) {
    row_index i = boost::python::extract<row_index>(ij[0]);
    column_index j = boost::python::extract<column_index>(ij[1]);
    self(i,j) = x;
  }

  static value_type getitem(wt& self, boost::python::tuple ij) {
    row_index i = boost::python::extract<row_index>(ij[0]);
    column_index j = boost::python::extract<column_index>(ij[1]);
    return self(i,j);
  }

  static void wrap() {
    using namespace boost::python;
    class_<wt>("matrix", no_init)
      .def(init<typename wt::row_index, typename wt::column_index>())
      .add_property("n_cols", &wt::n_cols)
      .add_property("n_rows", &wt::n_rows)
      .def("col",
           static_cast<vector<T>& (wt::*)(column_index)>(&wt::col),
           return_internal_reference<1>())
      .def("__getitem__", getitem)
      .def("__setitem__", setitem)
      .def("is_upper_triangular", &wt::is_upper_triangular)
      .def("is_unit_lower_triangular", &wt::is_unit_lower_triangular)
      .def("deep_copy", &wt::deep_copy)
      .def("transpose", &wt::transpose)
      .def(self*vector<T>())
    ;
  }

};

void wrap_matrix() {
  matrix_wrapper<double>::wrap();
}

}}}
