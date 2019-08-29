#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <scitbx/array_family/loops.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/boost_python/iterator_wrappers.h>

namespace scitbx { namespace af {

namespace {

  struct nested_loop_wrappers
  {
    typedef flex_grid_default_index_type array_t;
    typedef nested_loop<array_t> w_t;

    static array_t
    next(w_t& o)
    {
      if (o.over()) {
        PyErr_SetString(PyExc_StopIteration, "At end of loop.");
        boost::python::throw_error_already_set();
      }
      array_t result = o();
      o.incr();
      return result;
    }

    static void
    wrap()
    {
      using namespace boost::python;
      using boost::python::arg;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("nested_loop", no_init)
        .def(init<array_t const&, bool>(
          (arg("end"), arg("open_range")=true)))
        .def(init<array_t const&, array_t const&, bool>(
          (arg("begin"), arg("end"), arg("open_range")=true)))
        .def("incr", &w_t::incr)
        .def("begin", &w_t::begin, ccr())
        .def("end", &w_t::end, ccr())
        .def("__call__", &w_t::operator(), ccr())
        .def("over", &w_t::over)
        .def("__iter__", scitbx::boost_python::pass_through)
        .def("next", next)
        .def("__next__", next)
      ;
    }
  };

} // namespace <anoymous>

namespace boost_python {

  void wrap_loops()
  {
    nested_loop_wrappers::wrap();
  }

}}} // namespace scitbx::af::boost_python
