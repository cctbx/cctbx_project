#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/structure_factors_direct.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace xray { namespace structure_factors {
namespace boost_python {

namespace {

  struct direct_wrappers
  {
    typedef direct<> w_t;
    typedef w_t::scatterer_type scatterer_type;
    typedef w_t::float_type float_type;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("structure_factors_direct", no_init)
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<scatterer_type> const&>())
        .def(init<math::cos_sin_table<float_type> const&,
                  uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<scatterer_type> const&>())
        .def("f_calc", &w_t::f_calc)
      ;
    }
  };

} // namespace <anoymous>

}} // namespace structure_factors::boost_python

namespace boost_python {

  void wrap_structure_factors_direct()
  {
    structure_factors::boost_python::direct_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
