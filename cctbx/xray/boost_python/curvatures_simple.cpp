#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/curvatures_simple.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx { namespace xray { namespace structure_factors {
namespace curvatures_simple {
namespace boost_python {

namespace {

  struct grads_and_curvs_target_wrappers
  {
    typedef grads_and_curvs_target<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>(
          "structure_factors_curvatures_simple_grads_and_curvs_target",
          no_init)
        .def(init<
          uctbx::unit_cell const&,
          sgtbx::space_group const&,
          af::const_ref<scatterer<> > const&,
          xray::scattering_type_registry const&,
          sgtbx::site_symmetry_table const&,
          af::const_ref<miller::index<> > const&,
          af::const_ref<std::complex<double> > const&,
          af::const_ref<scitbx::vec3<double> > const&>((
            arg("unit_cell"),
            arg("space_group"),
            arg("scatterers"),
            arg("scattering_type_registry"),
            arg("site_symmetry_table"),
            arg("miller_indices"),
            arg("da_db"),
            arg("daa_dbb_dab"))))
        .add_property("grads", make_getter(&w_t::grads, rbv()))
        .add_property("curvs", make_getter(&w_t::curvs, rbv()))
      ;
    }
  };

} // namespace <anoymous>

}}} // namespace structure_factors::curvatures_simple::boost_python

namespace boost_python {

  void wrap_curvatures_simple()
  {
    structure_factors::curvatures_simple::boost_python
      ::grads_and_curvs_target_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
