#include <boost/python/class.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/str.hpp>
#include <boost/python/extract.hpp>

#include <scitbx/sparse/matrix.h>
#include <scitbx/sparse/io.h>

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

  static boost::python::str as_mathematica(wt const &m) {
    std::stringstream o(std::ios_base::out);
    o << dense_display(m);
    return boost::python::str(o.str().c_str());
  }

  static wt & permute_rows(wt &m, af::const_ref<row_index> const &p) {
    return m.permute_rows(p);
  }

  static void wrap() {
    using namespace boost::python;
    return_internal_reference<> rir;
    class_<wt>("matrix", no_init)
      .def(init<boost::optional<typename wt::row_index>,
                typename wt::column_index>())
      .add_property("n_cols", &wt::n_cols)
      .add_property("n_rows", &wt::n_rows)
      .def("col",
           static_cast<vector<T>& (wt::*)(column_index)>(&wt::col),
           return_internal_reference<1>())
      .def("__getitem__", getitem)
      .def("__setitem__", setitem)
      .def("is_structural_zero", &wt::is_structural_zero)
      .def("compact", &wt::compact)
      .def("is_upper_triangular", &wt::is_upper_triangular)
      .def("is_unit_lower_triangular", &wt::is_unit_lower_triangular)
      .def("deep_copy", &wt::deep_copy)
      .def("transpose", &wt::transpose)
      .def("permute_rows", permute_rows, arg("permutation"), rir)
      .def(self*vector<T>())
      .def(self*self)
        .def("self_transpose_times_symmetric_times_self",
           &wt::this_transpose_times_symmetric_times_this)
      .def(typename wt::dense_vector_const_ref() * self)
      .def("as_mathematica", as_mathematica)
    ;
  }

};

void wrap_matrix() {
  matrix_wrapper<double>::wrap();
}

}}}
