#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/structure_factors_simple.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace xray { namespace structure_factors {
namespace boost_python {

namespace {

  struct simple_wrappers
  {
    typedef simple<> w_t;
    typedef w_t::scatterer_type scatterer_type;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("structure_factors_simple", no_init)
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<scatterer_type> const&,
                  scattering_dictionary const&>())
        .def("f_calc", &w_t::f_calc)
      ;
    }
  };

} // namespace <anoymous>

}} // namespace structure_factors::boost_python

namespace boost_python {

  void wrap_structure_factors_simple()
  {
    structure_factors::boost_python::simple_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
