#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <cctbx/crystal/symmetry.h>

namespace cctbx { namespace crystal {

namespace {

  struct symmetry_wrappers : boost::python::pickle_suite
  {
    typedef symmetry w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.unit_cell(),
        self.space_group()
        );
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      class_<w_t>("symmetry", no_init)
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&>(
          (arg("unit_cell"), arg("space_group"))))
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("space_group", &w_t::space_group, rir())
        .def("change_basis", &w_t::change_basis,
          (arg("change_of_basis_op")))
        .def_pickle(symmetry_wrappers())
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
