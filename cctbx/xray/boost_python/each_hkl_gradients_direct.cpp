#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/each_hkl_gradients_direct.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace xray { namespace structure_factors {
namespace boost_python {

namespace {

  struct each_hkl_gradients_direct_wrappers
  {
    typedef each_hkl_gradients_direct<> w_t;
    typedef w_t::scatterer_type scatterer_type;
    typedef w_t::float_type float_type;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("each_hkl_gradients_direct", no_init)
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<scatterer_type> const&,
                  af::const_ref<float_type> const&,
                  scattering_type_registry const&,
                  sgtbx::site_symmetry_table const&,
                  std::size_t>())
        .def(init<math::cos_sin_table<double> const&,
                  uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<scatterer_type> const&,
                  af::const_ref<float_type> const&,
                  scattering_type_registry const&,
                  sgtbx::site_symmetry_table const&,
                  std::size_t>())
        .def("d_fcalc_d_fp", &w_t::d_fcalc_d_fp)
        .def("d_fcalc_d_fdp", &w_t::d_fcalc_d_fdp)
      ;
    }
  };

} // namespace <anoymous>

}} // namespace structure_factors::boost_python

namespace boost_python {

  void wrap_each_hkl_gradients_direct()
  {
    structure_factors::boost_python::each_hkl_gradients_direct_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
