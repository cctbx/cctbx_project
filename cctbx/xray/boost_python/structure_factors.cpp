#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/structure_factors.h>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>

namespace cctbx { namespace xray { namespace structure_factors {
namespace boost_python {

namespace {

  struct direct_with_first_derivatives_wrappers
  {
    typedef direct_with_first_derivatives<> w_t;
    typedef w_t::scatterer_type scatterer_type;
    typedef w_t::float_type float_type;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("structure_factors_direct_with_first_derivatives", no_init)
        .def(init<uctbx::unit_cell const&,
                  sgtbx::space_group const&,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<scatterer_type> const&,
                  af::const_ref<std::complex<float_type> > const&,
                  gradient_flags const&>())
        .def("f_calc", &w_t::f_calc)
        .def("d_target_d_site", &w_t::d_target_d_site)
        .def("d_target_d_u_iso", &w_t::d_target_d_u_iso)
        .def("d_target_d_u_star", &w_t::d_target_d_u_star)
        .def("d_target_d_occupancy", &w_t::d_target_d_occupancy)
        .def("d_target_d_fp", &w_t::d_target_d_fp)
        .def("d_target_d_fdp", &w_t::d_target_d_fdp)
      ;
    }
  };

  void
  wrap_functions()
  {
    using namespace boost::python;

    def("structure_factors_d_target_d_site_in_place_frac_as_cart",
      (void(*)(uctbx::unit_cell const&,
               af::ref<scitbx::vec3<double> > const&))
        d_target_d_site_in_place_frac_as_cart);
  }

} // namespace <anoymous>

}} // namespace structure_factors::boost_python

namespace boost_python {

  void wrap_structure_factors()
  {
    structure_factors::boost_python
    ::direct_with_first_derivatives_wrappers::wrap();

    structure_factors::boost_python::wrap_functions();
  }

}}} // namespace cctbx::xray::boost_python
