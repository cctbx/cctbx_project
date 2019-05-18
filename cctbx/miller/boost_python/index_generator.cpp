#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/index_generator.h>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <boost/python/class.hpp>
#include <boost/python/return_internal_reference.hpp>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(cctbx::miller::index_generator)

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct index_generator_wrappers
  {
    typedef index_generator w_t;

    static index<>
    next(w_t& o)
    {
      index<> result = o.next();
      if (result.is_zero()) {
        PyErr_SetString(PyExc_StopIteration, "At end of iteration.");
        boost::python::throw_error_already_set();
      }
      return result;
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      class_<w_t>("index_generator", no_init)
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group_type const&,
                  bool,
                  double>(
             args("unit_cell", "space_group_type", "anomalous_flag",
                  "resolution_limit")))
        .def(init<sgtbx::space_group_type const&,
                  bool,
                  index<> const&>(
             args("space_group_type", "anomalous_flag", "max_index")))
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("space_group_type", &w_t::space_group_type, rir())
        .def("anomalous_flag", &w_t::anomalous_flag)
        .def("asu", &w_t::asu, rir())
        .def("next", next)
        .def("__next__", next)
        .def("__iter__", scitbx::boost_python::pass_through)
        .def("to_array", &w_t::to_array)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_index_generator()
  {
    index_generator_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
