#include <boost/python/module.hpp>

namespace cctbx { namespace xray { namespace boost_python {

  void wrap_fast_gradients();
  void wrap_gradient_flags();
  void wrap_gradients_direct();
  void wrap_conversions();
  void wrap_minimization();
  void wrap_scatterer_flags();
  void wrap_sampling_base();
  void wrap_sampled_model_density();
  void wrap_scatterer();
  void wrap_scattering_type_registry();
  void wrap_structure_factors_direct();
  void wrap_structure_factors_raw_multithreaded_direct();
  void wrap_structure_factors_simple();
  void wrap_curvatures_simple();
  void wrap_targets();
  void wrap_f_model_core_data();
  void wrap_twin_targets();
  void wrap_grouped_data();
  void wrap_parameter_map();
  void wrap_twin_component();
  void wrap_extinction_correction();

namespace {

  void init_module()
  {
    wrap_conversions();
    wrap_gradient_flags();
    wrap_gradients_direct();
    wrap_sampling_base();
    wrap_fast_gradients();
    wrap_minimization();
    wrap_scatterer_flags();
    wrap_sampled_model_density();
    wrap_scatterer();
    wrap_scattering_type_registry();
    wrap_structure_factors_direct();
    #ifndef BOOST_DISABLE_THREADS
    wrap_structure_factors_raw_multithreaded_direct();
    #endif
    wrap_structure_factors_simple();
    wrap_curvatures_simple();
    wrap_targets();
    wrap_f_model_core_data();
    wrap_twin_targets();
    wrap_grouped_data();
    wrap_parameter_map();
    wrap_twin_component();
    wrap_extinction_correction();
  }

} // namespace <anonymous>
}}} // namespace cctbx::xray::boost_python

BOOST_PYTHON_MODULE(cctbx_xray_ext)
{
  cctbx::xray::boost_python::init_module();
}
