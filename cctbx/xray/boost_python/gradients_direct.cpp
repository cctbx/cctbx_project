#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/gradients_direct.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace xray { namespace structure_factors {
namespace boost_python {

namespace {

  struct gradients_direct_wrappers
  {
    typedef gradients_direct<> w_t;
    typedef w_t::scatterer_type scatterer_type;
    typedef w_t::float_type float_type;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("structure_factors_gradients_direct", no_init)
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<scatterer_type> const&,
                  af::const_ref<std::complex<float_type> > const&,
                  gradient_flags const&>())
        .def(init<math::cos_sin_table<double> const&,
                  uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<scatterer_type> const&,
                  af::const_ref<std::complex<float_type> > const&,
                  gradient_flags const&>())
        .def("d_target_d_site", &w_t::d_target_d_site)
        .def("d_target_d_u_iso", &w_t::d_target_d_u_iso)
        .def("d_target_d_u_star", &w_t::d_target_d_u_star)
        .def("d_target_d_occupancy", &w_t::d_target_d_occupancy)
        .def("d_target_d_fp", &w_t::d_target_d_fp)
        .def("d_target_d_fdp", &w_t::d_target_d_fdp)
      ;
    }
  };

} // namespace <anoymous>

}} // namespace structure_factors::boost_python

namespace boost_python {

  void wrap_gradients_direct()
  {
    structure_factors::boost_python::gradients_direct_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
