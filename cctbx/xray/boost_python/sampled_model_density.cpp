#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/sampled_model_density.h>
#include <boost/python/class.hpp>

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
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t, bases<w_t::base_t> >("sampled_model_density", no_init)
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<scatterer<> > const&,
                  scattering_dictionary const&,
                  grid_point_type const&,
                  grid_point_type const&,
                  optional<double const&,
                           double const&,
                           double const&,
                           bool,
                           bool,
                           double const&> >(
          (arg_("unit_cell"),
           arg_("scatterers"),
           arg_("scattering_dict"),
           arg_("fft_n_real"),
           arg_("fft_m_real"),
           arg_("u_base"),
           arg_("wing_cutoff"),
           arg_("exp_table_one_over_step_size"),
           arg_("force_complex"),
           arg_("sampled_density_must_be_positive"),
           arg_("tolerance_positive_definite"))))
        .def("real_map", &w_t::real_map)
        .def("complex_map", &w_t::complex_map)
        .def("eliminate_u_extra_and_normalize",
          &w_t::eliminate_u_extra_and_normalize,
          (arg_("miller_indices"), arg_("structure_factors")))
      ;
    }
  };

} // namespace <anoymous>

  void wrap_sampled_model_density()
  {
    sampled_model_density_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
