#include <boost/python/class.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/str.hpp>
#include <boost/python/slice.hpp>
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

  static boost::python::object setitem(wt& self,
                                       boost::python::tuple ij,
                                       boost::python::object x)
  {
    using namespace boost::python;
    object none;
    extract<row_index> p_i(ij[0]);
    extract<column_index> p_j(ij[1]);
    if (p_j.check()) {
      if (p_i.check()) {
        self(p_i(), p_j()) = extract<T>(x);
        return x;
      }
      else {
        extract<slice> p_i(ij[0]);
        if (p_i.check()) {
          slice i = p_i();
          if (i.start() == none && i.stop() == none) {
            self.col(p_j()) = extract< vector<T> >(x);
            return x;
          }
        }
      }
    }
                PyErr_SetString(PyExc_RuntimeError,
                    "Only self[i,j] = float() and self[:,j] = sparse.vector()"
                    " are supported.");
  }

  static boost::python::object getitem(wt& self, boost::python::tuple ij) {
    using namespace boost::python;
    object none;
    extract<row_index> p_i(ij[0]);
    extract<column_index> p_j(ij[1]);
    if (p_j.check()) {
      if (p_i.check()) {
        T x = self(p_i(), p_j());
        return object(x);
      }
      else {
        extract<slice> p_i(ij[0]);
        if (p_i.check()) {
          slice i = p_i();
          if (i.start() == none && i.stop() == none) {
            return object(self.col(p_j()));
          }
        }
      }
    }
                PyErr_SetString(PyExc_RuntimeError,
                    "Only self[i,j] and self[:,j] are supported.");
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
      .def(init<typename wt::row_index,
                typename wt::column_index>
           ((arg("rows"), arg("columns"))))
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
      .def(self*vector<T>(0))
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
