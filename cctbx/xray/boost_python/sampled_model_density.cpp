#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/sampled_model_density.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct sampled_model_density_wrappers
  {
    typedef sampled_model_density<> w_t;
    typedef w_t::grid_point_type grid_point_type;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t, bases<w_t::base_t> >("sampled_model_density", no_init)
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<scatterer<> > const&,
                  scattering_type_registry const&,
                  grid_point_type const&,
                  grid_point_type const&,
                  optional<double const&,
                           double const&,
                           double const&,
                           bool,
                           bool,
                           double const&,
                           bool,
                           int> >(
          (arg("unit_cell"),
           arg("scatterers"),
           arg("scattering_type_registry"),
           arg("fft_n_real"),
           arg("fft_m_real"),
           arg("u_base")=0.25,
           arg("wing_cutoff")=1.e-3,
           arg("exp_table_one_over_step_size")=-100,
           arg("force_complex")=false,
           arg("sampled_density_must_be_positive")=false,
           arg("tolerance_positive_definite")=1.e-5,
           arg("use_u_base_as_u_extra")=false,
           arg("store_grid_indices_for_each_scatterer")=0)))
        .def("real_map", &w_t::real_map)
        .def("complex_map", &w_t::complex_map)
        .def("eliminate_u_extra_and_normalize",
          &w_t::eliminate_u_extra_and_normalize,
          (arg("miller_indices"), arg("structure_factors")))
        .def("grid_indices_for_each_scatterer",
          &w_t::grid_indices_for_each_scatterer, ccr())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_sampled_model_density()
  {
    sampled_model_density_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
