#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/scatterer.h>
#include <cctbx/xray/scatterer_utils.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
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
    apply_symmetry_overloads, apply_symmetry, 3, 7)

  // work around Visual C++ 7 internal compiler error
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    structure_factor_array_overloads, structure_factor_array, 4, 4)

namespace {

  void init_module()
  {
    using namespace boost::python;

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

    def("apply_symmetry",
      (af::shared<std::size_t>(*)(
        uctbx::unit_cell const&,
        sgtbx::space_group const&,
        af::ref<scatterer<> > const&,
        double, double, bool, bool)) 0, apply_symmetry_overloads());

    def("rotate",
      (af::shared<scatterer<> >(*)(
        uctbx::unit_cell const&,
        scitbx::mat3<double> const&,
        af::const_ref<scatterer<> > const&)) rotate);
  }

} // namespace <anonymous>
}}} // namespace cctbx::xray::boost_python

BOOST_PYTHON_MODULE(xray_ext)
{
  cctbx::xray::boost_python::init_module();
}
