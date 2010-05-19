#include <boost/python/class.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/str.hpp>
#include <boost/python/slice.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/make_constructor.hpp>

#include <scitbx/sparse/matrix.h>
#include <scitbx/sparse/io.h>
#include <scitbx/sparse/boost_python/vector.h>

namespace scitbx { namespace sparse { namespace boost_python {

template<typename T>
struct matrix_wrapper
{
  typedef matrix<T> wt;
  typedef typename wt::index_type index_type;
  typedef typename wt::value_type value_type;

  static wt *from_list_of_dict(index_type m, index_type n,
                               boost::python::list cols)
  {
    using namespace boost::python;
    SCITBX_ASSERT(len(cols) == n);
    wt *result = new wt(m, n);
    for (index_type j=0; j<n; ++j) {
      dict col = extract<dict>(cols[j]);
      result->col(j) = vector_from_dict<T>(m, col);
    }
    return result;
  }

  static boost::python::object setitem(wt& self,
                                       boost::python::tuple ij,
                                       boost::python::object x)
  {
    using namespace boost::python;
    object none;
    extract<index_type> p_i(ij[0]);
    extract<index_type> p_j(ij[1]);
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
    else {
      throw scitbx::error("Only self[i,j] = float() "
                          "and self[:,j] = sparse.vector() are supported.");
    }
  }

  static boost::python::object getitem(wt& self, boost::python::tuple ij) {
    using namespace boost::python;
    object none;
    extract<index_type> p_i(ij[0]);
    extract<index_type> p_j(ij[1]);
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
    else {
      throw scitbx::error("Only self[i,j] and self[:,j] are supported.");
    }
  }

  static boost::python::str str_(wt const &m) {
    std::stringstream o(std::ios_base::out);
    o << dense_display(m);
    return boost::python::str(o.str().c_str());
  }

  static boost::python::str repr(wt const &m) {
    std::stringstream o(std::ios_base::out);
    std::string start("sparse.matrix(");
    o << start << "rows=" << m.n_rows() << ", columns=" << m.n_cols() << ",\n";
    std::string elts("elements_by_columns=[ ");
    o << std::setw(start.size()) << "" << elts;
    for (index_type j=0; j<m.n_cols(); ++j) {
      if (j > 0) o << std::setw(start.size() + elts.size()) << "";
      o << compressed_display(m.col(j)) << ",";
      if (j+1 < m.n_cols()) o << "\n";
    }
    o << " ])";
    return boost::python::str(o.str().c_str());
  }

  static wt & permute_rows(wt &m, af::const_ref<index_type> const &p) {
    return m.permute_rows(p);
  }

  static void wrap() {
    using namespace boost::python;
    return_internal_reference<> rir;
    class_<wt>("matrix", no_init)
      .def(init<typename wt::index_type,
                typename wt::index_type>
           ((arg("rows"), arg("columns"))))
      .def("__init__",
           make_constructor(from_list_of_dict,
                            default_call_policies(),
                            (arg("rows"), arg("columns"),
                             arg("elements_by_columns"))))
      .add_property("n_cols", &wt::n_cols)
      .add_property("n_rows", &wt::n_rows)
      .def("col",
           static_cast<vector<T>& (wt::*)(index_type)>(&wt::col),
           return_internal_reference<1>())
      .def("__getitem__", getitem)
      .def("__setitem__", setitem)
      .def("is_structural_zero", &wt::is_structural_zero)
      .def("compact", &wt::compact)
      .def("is_upper_triangular", &wt::is_upper_triangular)
      .def("is_unit_lower_triangular", &wt::is_unit_lower_triangular)
      .add_property("non_zeroes", &wt::non_zeroes)
      .def("deep_copy", &wt::deep_copy)
      .def("transpose", &wt::transpose)
      .def("permute_rows", permute_rows, arg("permutation"), rir)
      .def(self*vector<T>(0))
      .def(self * typename wt::dense_vector_const_ref())
      .def(typename wt::dense_vector_const_ref() * self)
      .def(self*self)
      .def("self_transpose_times_symmetric_times_self",
           &wt::this_transpose_times_symmetric_times_this)
      .def("__str__", str_)
      .def("__repr__", repr)
    ;
  }

};

void wrap_matrix() {
  matrix_wrapper<double>::wrap();
}

}}}
