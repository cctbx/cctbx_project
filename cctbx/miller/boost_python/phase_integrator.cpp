#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/phase_integrator.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct phase_integrator_wrappers
  {
    typedef phase_integrator w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("phase_integrator", no_init)
        .def(init<optional<unsigned> >(arg_("n_steps")=360/5))
        .def("n_steps", &w_t::n_steps)
        .def("__call__",
          (std::complex<double>(w_t::*)(
            sgtbx::phase_info const&,
            hendrickson_lattman<> const&) const) &w_t::operator(), (
          arg_("phase_info"), arg_("hendrickson_lattman")))
        .def("__call__",
          (af::shared<std::complex<double> >(w_t::*)(
            sgtbx::space_group const&,
            af::const_ref<miller::index<> > const&,
            af::const_ref<hendrickson_lattman<> > const&) const)
              &w_t::operator(), (
          arg_("space_group"),
          arg_("miller_indices"),
          arg_("hendrickson_lattman_coefficients")))
      ;
    }
  };

} // namespace <anoymous>

  void wrap_phase_integrator()
  {
    phase_integrator_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
