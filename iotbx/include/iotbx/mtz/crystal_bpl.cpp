#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <iotbx/mtz/dataset.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

namespace iotbx { namespace mtz {
namespace {

  struct crystal_wrappers
  {
    typedef crystal w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("crystal", no_init)
        .def(init<object const&, int>((arg_("mtz_object"), arg_("i_crystal"))))
        .def("mtz_object", &w_t::mtz_object)
        .def("i_crystal", &w_t::i_crystal)
        .def("id", &w_t::id)
        .def("name", &w_t::name)
        .def("project_name", &w_t::project_name)
        .def("unit_cell_parameters", &w_t::unit_cell_parameters)
        .def("unit_cell", &w_t::unit_cell)
        .def("n_datasets", &w_t::n_datasets)
        .def("datasets", &w_t::datasets)
        .def("add_dataset", &w_t::add_dataset, (
          arg_("name"), arg_("wavelength")))
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_crystal");
      }
    }
  };

  void
  wrap_all()
  {
    crystal_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_crystal() { wrap_all(); }

}}} // namespace iotbx::mtz::boost_python
