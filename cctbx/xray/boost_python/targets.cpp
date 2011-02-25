#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/targets.h>
#include <cctbx/xray/targets/least_squares.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace xray { namespace targets { namespace boost_python {

namespace {

  struct common_results_wrappers
  {
    typedef common_results w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("targets_common_results", no_init)
        .def(init<
          af::shared<double> const&,
          double,
          boost::optional<double> const&,
          af::shared<std::complex<double> > const&>((
            arg("target_per_reflection"),
            arg("target_work"),
            arg("target_test"),
            arg("gradients_work"))))
        .def(init<
          af::shared<double> const&,
          double,
          boost::optional<double> const&,
          af::shared<std::complex<double> > const&,
          af::shared<scitbx::vec3<double> > const&>((
            arg("target_per_reflection"),
            arg("target_work"),
            arg("target_test"),
            arg("gradients_work"),
            arg("hessians_work"))))
        .def("target_per_reflection", &w_t::target_per_reflection, ccr())
        .def("target_work", &w_t::target_work)
        .def("target_test", &w_t::target_test)
        .def("gradients_work", &w_t::gradients_work, ccr())
        .def("hessians_work", &w_t::hessians_work, ccr())
      ;
    }
  };

  struct least_squares_wrappers
  {
    typedef least_squares w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<common_results> >(
          "targets_least_squares", no_init)
        .def(init<
          bool,
          char,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          af::const_ref<bool> const&,
          af::const_ref< std::complex<double> > const&,
          int,
          double>((
            arg("compute_scale_using_all_data"),
            arg("obs_type"),
            arg("obs"),
            arg("weights"),
            arg("r_free_flags"),
            arg("f_calc"),
            arg("derivatives_depth"),
            arg("scale_factor"))))
        .def("compute_scale_using_all_data",&w_t::compute_scale_using_all_data)
        .def("obs_type", &w_t::obs_type)
        .def("scale_factor", &w_t::scale_factor)
      ;
    }
  };

  template <template<typename> class FcalcFunctor>
  struct least_squares_residual_wrappers
  {
    typedef least_squares_residual<FcalcFunctor> w_t;

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      class_<w_t>(python_name,
                  "Boost.Python wrapping of the C++ class"
                  "U{least_squares_residual<c_plus_plus/"
                  "classcctbx_1_1xray_1_1targets_1_1least__squares__residual.html>}",
                  no_init)
        .def(init<af::const_ref<double> const&,
                  af::const_ref<double> const&,
                  af::const_ref<std::complex<double> > const&,
                  optional<bool, double> >())
        .def(init<af::const_ref<double> const&,
                  af::const_ref<std::complex<double> > const&,
                  optional<bool, double> >())
        .def("scale_factor", &w_t::scale_factor)
        .def("target", &w_t::target)
        .def("derivatives", &w_t::derivatives)
      ;
    }
  };

  struct intensity_correlation_wrappers
  {
    typedef intensity_correlation<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("targets_intensity_correlation", no_init)
        .def(init<af::const_ref<double> const&,
                  af::const_ref<int> const&,
                  af::const_ref<std::complex<double> > const&,
                  optional<bool> >())
        .def(init<af::const_ref<double> const&,
                  af::const_ref<std::complex<double> > const&,
                  optional<bool> >())
        .def("correlation", &w_t::correlation)
        .def("target", &w_t::target)
        .def("derivatives", &w_t::derivatives)
      ;
    }
  };

  struct maximum_likelihood_criterion_wrappers
  {
    typedef maximum_likelihood_criterion w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<common_results> >(
          "targets_maximum_likelihood_criterion", no_init)
        .def(init<
          af::const_ref<double> const&,
          af::const_ref<bool> const&,
          af::const_ref<std::complex<double> > const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          double,
          af::const_ref<int> const&,
          af::const_ref<bool> const&,
          bool>((
            arg("f_obs"),
            arg("r_free_flags"),
            arg("f_calc"),
            arg("alpha"),
            arg("beta"),
            arg("scale_factor"),
            arg("epsilons"),
            arg("centric_flags"),
            arg("compute_gradients"))))
      ;
    }
  };

  struct maximum_likelihood_criterion_hl_wrappers
  {
    typedef maximum_likelihood_criterion_hl w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<common_results> >(
          "targets_maximum_likelihood_criterion_hl", no_init)
        .def(init<
          af::const_ref<double> const&,
          af::const_ref<bool>,
          af::const_ref<cctbx::hendrickson_lattman<double> > const&,
          af::const_ref<std::complex<double> > const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          af::const_ref<int> const&,
          af::const_ref<bool> const&,
          double,
          bool>((
            arg("f_obs"),
            arg("r_free_flags"),
            arg("experimental_phases"),
            arg("f_calc"),
            arg("alpha"),
            arg("beta"),
            arg("epsilons"),
            arg("centric_flags"),
            arg("integration_step_size"),
            arg("compute_gradients"))))
      ;
    }
  };

} // namespace <anoymous>

}} // namespace targets::boost_python

namespace boost_python {

  void wrap_targets()
  {
    targets::boost_python::common_results_wrappers::wrap();
    targets::boost_python::least_squares_wrappers::wrap();
    targets::boost_python::least_squares_residual_wrappers<
      cctbx::xray::targets::f_calc_modulus>::wrap(
      "targets_least_squares_residual");
    targets::boost_python::least_squares_residual_wrappers<
      cctbx::xray::targets::f_calc_modulus_square>::wrap(
      "targets_least_squares_residual_for_intensity");
    targets::boost_python::intensity_correlation_wrappers::wrap();
    targets::boost_python::maximum_likelihood_criterion_wrappers::wrap();
    targets::boost_python::maximum_likelihood_criterion_hl_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
