#include <boost/python/class.hpp>
#include <cctbx/xray/gradient_flags.h>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct gradient_flags_wrappers
  {
    typedef gradient_flags w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("gradient_flags", no_init)
        .def(init<bool, bool, bool, bool, bool, bool>())
        .def_readwrite("site", &w_t::site)
        .def_readwrite("u_iso", &w_t::u_iso)
        .def_readwrite("u_aniso", &w_t::u_aniso)
        .def_readwrite("occupancy", &w_t::occupancy)
        .def_readwrite("fp", &w_t::fp)
        .def_readwrite("fdp", &w_t::fdp)
        .def("all_false", &w_t::all_false)
        .def("adjust", &w_t::adjust)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_gradient_flags()
  {
    gradient_flags_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
