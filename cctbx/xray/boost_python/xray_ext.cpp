#include <boost/python/module.hpp>

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

namespace {

  void init_module()
  {
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
  }

} // namespace <anonymous>
}}} // namespace cctbx::xray::boost_python

BOOST_PYTHON_MODULE(cctbx_xray_ext)
{
  cctbx::xray::boost_python::init_module();
}
