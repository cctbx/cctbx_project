#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/scatterer_utils.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

namespace cctbx { namespace xray { namespace boost_python {

  void wrap_fast_gradients();
  void wrap_gradient_flags();
  void wrap_gradients_direct();
  void wrap_conversions();
  void wrap_minimization();
  void wrap_sampling_base();
  void wrap_sampled_model_density();
  void wrap_scatterer();
  void wrap_scattering_dictionary();
  void wrap_structure_factors_direct();
  void wrap_structure_factors_simple();
  void wrap_targets();

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    apply_symmetry_u_star_overloads, apply_symmetry_u_star, 3, 6)

namespace {

  void init_module()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_; // gcc 2.96 workaround

    wrap_conversions();
    wrap_gradient_flags();
    wrap_gradients_direct();
    wrap_sampling_base();
    wrap_fast_gradients();
    wrap_minimization();
    wrap_sampled_model_density();
    wrap_scatterer();
    wrap_scattering_dictionary();
    wrap_structure_factors_direct();
    wrap_structure_factors_simple();
    wrap_targets();

    def("apply_symmetry_site",
      (void(*)(
        sgtbx::site_symmetry_table const&,
        af::ref<scatterer<> > const&)) apply_symmetry_site, (
          arg_("site_symmetry_table"),
          arg_("scatterers")));

    def("apply_symmetry_u_star",
      (void(*)(
        uctbx::unit_cell const&,
        sgtbx::site_symmetry_table const&,
        af::ref<scatterer<> > const&,
        double, bool, bool)) 0, apply_symmetry_u_star_overloads((
          arg_("unit_cell"),
          arg_("site_symmetry_table"),
          arg_("scatterers"),
          arg_("u_star_tolerance")=0,
          arg_("assert_is_positive_definite")=false,
          arg_("assert_min_distance_sym_equiv")=true)));

    def("add_scatterers_ext",
      (void(*)(
        uctbx::unit_cell const&,
        sgtbx::space_group const&,
        af::ref<scatterer<> > const&,
        sgtbx::site_symmetry_table&,
        sgtbx::site_symmetry_table const&,
        double, double, bool, bool)) add_scatterers_ext, (
          arg_("unit_cell"),
          arg_("space_group"),
          arg_("scatterers"),
          arg_("site_symmetry_table"),
          arg_("site_symmetry_table_for_new"),
          arg_("min_distance_sym_equiv"),
          arg_("u_star_tolerance"),
          arg_("assert_is_positive_definite"),
          arg_("assert_min_distance_sym_equiv")));

    def("rotate",
      (af::shared<scatterer<> >(*)(
        uctbx::unit_cell const&,
        scitbx::mat3<double> const&,
        af::const_ref<scatterer<> > const&)) rotate, (
          arg_("unit_cell"),
          arg_("rotation_matrix"),
          arg_("scatterers")));
  }

} // namespace <anonymous>
}}} // namespace cctbx::xray::boost_python

BOOST_PYTHON_MODULE(cctbx_xray_ext)
{
  cctbx::xray::boost_python::init_module();
}
