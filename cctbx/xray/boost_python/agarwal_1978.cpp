#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/agarwal_1978.h>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct agarwal_1978_wrappers
  {
    typedef agarwal_1978<> w_t;
    typedef w_t::grid_point_type grid_point_type;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("agarwal_1978", no_init)
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<scatterer<> > const&,
                  af::const_ref<std::complex<double>,
                                maptbx::c_grid_padded_p1<3> > const&,
                  bool,
                  grid_point_type const&,
                  grid_point_type const&,
                  optional<double const&,
                           double const&,
                           double const&,
                           bool,
                           bool> >())
        .def("unit_cell", &w_t::unit_cell, rir())
        .def("u_extra", &w_t::u_extra)
        .def("wing_cutoff", &w_t::wing_cutoff)
        .def("exp_table_one_over_step_size",
          &w_t::exp_table_one_over_step_size)
        .def("n_scatterers_passed", &w_t::n_scatterers_passed)
        .def("n_contributing_scatterers", &w_t::n_contributing_scatterers)
        .def("n_anomalous_scatterers", &w_t::n_anomalous_scatterers)
        .def("anomalous_flag", &w_t::anomalous_flag)
        .def("exp_table_size", &w_t::exp_table_size)
        .def("max_shell_radii", &w_t::max_shell_radii, ccr())
        .def("max_shell_radii_frac", &w_t::max_shell_radii_frac)
        .def("grad", &w_t::grad)
        .def("grad_x", &w_t::grad_x)
        .def("grad_y", &w_t::grad_y)
        .def("grad_z", &w_t::grad_z)
        .def("grad_u_iso", &w_t::grad_u_iso)
        .def("grad_u_00", &w_t::grad_u_00)
        .def("grad_u_11", &w_t::grad_u_11)
        .def("grad_u_22", &w_t::grad_u_22)
        .def("grad_u_01", &w_t::grad_u_01)
        .def("grad_u_02", &w_t::grad_u_02)
        .def("grad_u_12", &w_t::grad_u_12)
        .def("grad_occupancy", &w_t::grad_occupancy)
        .def("grad_fp", &w_t::grad_fp)
      ;
    }
  };

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    calc_u_extra_overloads, calc_u_extra, 2, 4)

} // namespace <anoymous>

  void wrap_agarwal_1978()
  {
    using namespace boost::python;
    agarwal_1978_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
