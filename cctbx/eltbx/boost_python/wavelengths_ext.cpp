#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <cctbx/eltbx/wavelengths.h>

namespace cctbx { namespace eltbx { namespace wavelengths {
namespace boost_python {

namespace {

  struct characteristic_wrappers
  {
    typedef characteristic w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("characteristic", no_init)
        .def(init<std::string const&>())
        .def("label", &w_t::label)
        .def("as_angstrom", &w_t::as_angstrom)
        .def("as_kev", &w_t::as_kev)
        .def("as_ev", &w_t::as_ev)
      ;
    }
  };

  void init_module()
  {
    characteristic_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      characteristic,
      characteristic_iterator>::wrap(
        "characteristic_iterator");
  }

} // namespace <anonymous>
}}}} // namespace cctbx::eltbx::wavelengths::boost_python

BOOST_PYTHON_MODULE(cctbx_eltbx_wavelengths_ext)
{
  cctbx::eltbx::wavelengths::boost_python::init_module();
}
