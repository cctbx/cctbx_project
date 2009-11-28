#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <cctbx/sgtbx/space_group.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct phase_info_wrappers
  {
    typedef phase_info w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("phase_info", no_init)
        .def(init<
          space_group const&,
          miller::index<> const&,
          optional<bool> >((
            arg("space_group"),
            arg("miller_index"),
            arg("no_test_sys_absent")=false)))
        .def("sys_abs_was_tested", &w_t::sys_abs_was_tested)
        .def("is_sys_absent", &w_t::is_sys_absent)
        .def("is_centric", &w_t::is_centric)
        .def("ht", &w_t::ht)
        .def("t_den", &w_t::t_den)
        .def("ht_angle", &w_t::ht_angle, (arg("deg")=false))
        .def("is_valid_phase", &w_t::is_valid_phase, (
          arg("phi"),
          arg("deg")=false,
          arg("tolerance")=1e-5))
        .def("nearest_valid_phase",
          (double(w_t::*)(double, bool) const) &w_t::nearest_valid_phase, (
            arg("phi"),
            arg("deg")=false))
        .def("valid_structure_factor",
          (std::complex<double>(w_t::*)(std::complex<double> const&) const)
            &w_t::valid_structure_factor, (arg("f")))
      ;
    }
  };

} // namespace <anoymous>

  void wrap_phase_info()
  {
    phase_info_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
