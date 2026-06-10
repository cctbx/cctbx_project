#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/python/import.hpp>
#include <boost/python/str.hpp>
#include <cctbx/uctbx/niggli_reduction.h>

namespace cctbx { namespace uctbx { namespace boost_python {

namespace {

  struct niggli_reduction_wrappers
  {
    typedef niggli_reduction<> w_t;

    // Return r_inv as scitbx.matrix.sqr directly, bypassing the intermediate
    // boost.python mat3<int>->tuple conversion and the Python inject_into layer.
    // scitbx.matrix is guaranteed to be in sys.modules by the time any
    // niggli_reduction object is constructed (cctbx.array_family.flex imports
    // it transitively before cctbx_uctbx_ext is initialised).
    static boost::python::object
    change_of_basis_op_as_py(w_t const& self)
    {
      namespace bp = boost::python;
      bp::import("cctbx.sgtbx");
      return bp::object(self.change_of_basis_op());
    }

    static boost::python::object
    r_inv_as_sqr(w_t const& self)
    {
      namespace bp = boost::python;
      scitbx::mat3<int> const& m = self.r_inv();
      return bp::import("scitbx.matrix").attr("sqr")(
        bp::make_tuple(m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8]));
    }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("niggli_reduction", no_init)
        .def(init<unit_cell const&, double, std::size_t>((
          arg("unit_cell"),
          arg("relative_epsilon"),
          arg("iteration_limit"))))
        .def("as_gruber_matrix",         &w_t::as_gruber_matrix)
        .def("as_niggli_matrix",         &w_t::as_niggli_matrix)
        .def("as_sym_mat3",              &w_t::as_sym_mat3)
        .def("as_unit_cell",             &w_t::as_unit_cell)
        .def("change_of_basis_op",       change_of_basis_op_as_py)
        .def("r_inv",                    r_inv_as_sqr)
        .def("n_iterations",             &w_t::n_iterations)
        .def("iteration_limit",          &w_t::iteration_limit)
        .def("type",                     &w_t::type)
        .def("meets_primary_conditions", &w_t::meets_primary_conditions)
        .def("meets_main_conditions",    &w_t::meets_main_conditions)
        .def("is_buerger_cell",          &w_t::is_buerger_cell)
        .def("is_niggli_cell",           &w_t::is_niggli_cell)
      ;
    }
  };

  // Exception translator: maps niggli_reduction_iteration_limit_exceeded to
  // the Python class reduction_base.iteration_limit_exceeded.  This is
  // required so that client code using "except iteration_limit_exceeded:"
  // (imported from cctbx.uctbx.reduction_base) correctly catches exceptions
  // raised by the C++ algorithm.
  //
  // The import is deferred to raise-time (not module-init time) to avoid
  // circular imports: cctbx_uctbx_ext is imported by cctbx.uctbx.__init__,
  // which imports cctbx.uctbx.reduction_base, so the latter is guaranteed
  // to be loaded before any Niggli reduction can be called.
  void translate_niggli_iter_limit(
    niggli_reduction_iteration_limit_exceeded const& e)
  {
    namespace bp = boost::python;
    try {
      bp::object rb = bp::import("cctbx.uctbx.reduction_base");
      bp::object exc_class = rb.attr("iteration_limit_exceeded");
      PyErr_SetObject(exc_class.ptr(), bp::str(e.what()).ptr());
    } catch (...) {
      // Fallback: raise as a plain RuntimeError so something reaches Python.
      PyErr_SetString(PyExc_RuntimeError, e.what());
    }
  }

} // namespace <anonymous>

  void wrap_niggli_reduction()
  {
    boost::python::register_exception_translator<
      niggli_reduction_iteration_limit_exceeded>(
        translate_niggli_iter_limit);
    niggli_reduction_wrappers::wrap();
  }

}}} // namespace cctbx::uctbx::boost_python
