#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <cctbx/crystal/symmetry.h>

namespace cctbx { namespace crystal {

namespace {

  struct symmetry_wrappers
  {
    typedef symmetry w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_internal_reference<> rir;
      class_<w_t>("symmetry", no_init)
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&>(
          (arg_("unit_cell"), arg_("space_group"))))
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("space_group", &w_t::space_group, rir())
        .def("change_basis", &w_t::change_basis,
          (arg_("change_of_basis_op")))
      ;
    }
  };

} // namespace <anoymous>

namespace boost_python {

  void wrap_symmetry()
  {
    symmetry_wrappers::wrap();
  }

}}} // namespace cctbx::crystal::boost_python
