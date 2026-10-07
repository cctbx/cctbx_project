#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/targets.h>
#include <cctbx/xray/targets/least_squares.h>
#include <cctbx/xray/targets/correlation.h>
#include <cctbx/xray/targets/mlf.h>
#include <cctbx/xray/targets/mli.h>
#include <cctbx/xray/targets/mlhl.h>
#include <cctbx/xray/targets/llgi.h>
#include <cctbx/xray/targets/llgi_e.h>
#include <cctbx/xray/targets/llgi_exact.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_by_value.hpp>

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
          boost::optional<double> const&,
          boost::optional<double> const&,
          af::shared<std::complex<double> > const&>((
            arg("target_per_reflection"),
            arg("target_work"),
            arg("target_test"),
            arg("gradients_work"))))
        .def(init<
          af::shared<double> const&,
          boost::optional<double> const&,
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
        .def("target", &w_t::target_work) // backward compatibility
        .def("target_test", &w_t::target_test)
        .def("gradients_work", &w_t::gradients_work, ccr())
        .def("derivatives", &w_t::gradients_work, ccr())//backward compatibility
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
          af::const_ref<std::complex<double> > const&,
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

  struct correlation_wrappers
  {
    typedef correlation w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<common_results> >(
          "targets_correlation", no_init)
        .def(init<
          char,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          af::const_ref<bool> const&,
          af::const_ref<std::complex<double> > const&,
          int>((
            arg("obs_type"),
            arg("obs"),
            arg("weights"),
            arg("r_free_flags"),
            arg("f_calc"),
            arg("derivatives_depth"))))
        .def("obs_type", &w_t::obs_type)
        .def("cc", &w_t::cc)
        .def("correlation", &w_t::cc)
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

  struct mlf_wrappers
  {
    typedef mlf::target_and_gradients w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<common_results> >(
          "mlf_target_and_gradients", no_init)
        .def(init<
          af::const_ref<double> const&,
          af::const_ref<bool> const&,
          af::const_ref<std::complex<double> > const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          double,
          af::const_ref<double> const&,
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

  struct mli_wrappers
  {
    typedef mli::target_and_gradients w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<common_results> >(
          "mli_target_and_gradients", no_init)
        .def(init<
          af::const_ref<double> const&,
          af::const_ref<bool> const&,
          af::const_ref<std::complex<double> > const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          double,
          af::const_ref<double> const&,
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

  struct mlhl_wrappers
  {
    typedef mlhl::target_and_gradients w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<common_results> >(
          "mlhl_target_and_gradients", no_init)
        .def(init<
          af::const_ref<double> const&,
          af::const_ref<bool>,
          af::const_ref<cctbx::hendrickson_lattman<double> > const&,
          af::const_ref<std::complex<double> > const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
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

  struct llgi_wrappers
  {
    typedef llgi::target_and_gradients w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<common_results> >(
          "llgi_target_and_gradients", no_init)
        .def(init<
          af::const_ref<double> const&,
          af::const_ref<bool> const&,
          af::const_ref<std::complex<double> > const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          double,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          af::const_ref<bool> const&,
          bool,
          optional<llgi_exact::hybrid const*> >((
            arg("f_eff"),
            arg("r_free_flags"),
            arg("f_calc"),
            arg("dobs"),
            arg("sigmaa"),
            arg("scatfrac"),
            arg("scale_factor"),
            arg("teps"),
            arg("resn"),
            arg("centric_flags"),
            arg("compute_gradients"),
            arg("hybrid")=object())))
        .def("n_exact", &w_t::n_exact)
      ;
    }
  };

  struct llgi_sigmaa_scatfrac_wrappers
  {
    typedef llgi::sigmaa_scatfrac_target_and_gradients w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>(
          "llgi_sigmaa_scatfrac_target_and_gradients", no_init)
        .def(init<
          af::const_ref<double> const&,
          af::const_ref<bool> const&,
          af::const_ref<std::complex<double> > const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          double,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          af::const_ref<bool> const&,
          optional<llgi_exact::hybrid const*> >((
            arg("f_eff"),
            arg("selection"),
            arg("f_calc"),
            arg("dobs"),
            arg("sigmaa"),
            arg("scatfrac"),
            arg("scale_factor"),
            arg("teps"),
            arg("resn"),
            arg("centric_flags"),
            arg("hybrid")=object())))
        .def("target", &w_t::target)
        .def("n_exact", &w_t::n_exact)
        .def("d_target_by_dsigmaa", &w_t::d_target_by_dsigmaa, ccr())
        .def("d_target_by_dscatfrac", &w_t::d_target_by_dscatfrac, ccr())
      ;
    }
  };

  struct llgi_e_sigmaa_wrappers
  {
    typedef llgi_e::sigmaa_target_and_gradients w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>(
          "llgi_e_sigmaa_target_and_gradients", no_init)
        .def(init<
          af::const_ref<double> const&,
          af::const_ref<bool> const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          af::const_ref<bool> const&,
          optional<llgi_exact::hybrid const*> >((
            arg("e_eff"),
            arg("selection"),
            arg("e_model"),
            arg("dobs"),
            arg("sigmaa"),
            arg("centric_flags"),
            arg("hybrid")=object())))
        .def("target", &w_t::target)
        .def("n_exact", &w_t::n_exact)
        .def("d_target_by_dsigmaa", &w_t::d_target_by_dsigmaa, ccr())
      ;
    }
  };

  struct llgi_exact_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      {
        typedef llgi_exact::evaluate_many w_t;
        class_<w_t>("llgi_exact_evaluate", no_init)
          .def(init<
            af::const_ref<double> const&,
            af::const_ref<double> const&,
            af::const_ref<double> const&,
            af::const_ref<double> const&,
            af::const_ref<bool> const&>((
              arg("e_obs_sq"),
              arg("sig_e_obs_sq"),
              arg("e_calc"),
              arg("sigmaa"),
              arg("centric_flags"))))
          .def(init<
            af::const_ref<double> const&,
            af::const_ref<double> const&,
            af::const_ref<double> const&,
            af::const_ref<double> const&,
            af::const_ref<bool> const&,
            af::const_ref<double> const&>((
              arg("e_obs_sq"),
              arg("sig_e_obs_sq"),
              arg("e_calc"),
              arg("sigmaa"),
              arg("centric_flags"),
              arg("null_log_z"))))
          .add_property("ll", make_getter(&w_t::ll, return_value_policy<return_by_value>()))
          .add_property("d_ll_d_ec", make_getter(&w_t::d_ll_d_ec, return_value_policy<return_by_value>()))
          .add_property("d_ll_d_a", make_getter(&w_t::d_ll_d_a, return_value_policy<return_by_value>()))
          .add_property("e_expected", make_getter(&w_t::e_expected, return_value_policy<return_by_value>()))
          .add_property("e_abs_expected", make_getter(&w_t::e_abs_expected, return_value_policy<return_by_value>()))
        ;
      }
      {
        typedef llgi_exact::hybrid w_t;
        class_<w_t>("llgi_hybrid", no_init)
          .def(init<
            af::const_ref<double> const&,
            af::const_ref<double> const&,
            af::const_ref<bool> const&,
            af::const_ref<bool> const&,
            double>((
              arg("e_obs_sq"),
              arg("sig_e_obs_sq"),
              arg("force_exact"),
              arg("centric_flags"),
              arg("rice_kappa"))))
          .def("size", &w_t::size)
          .def("select", &w_t::select, (arg("selection")))
          .def("with_rice_kappa", &w_t::with_rice_kappa, (arg("t")))
          .def("exact_selection", &w_t::exact_selection, (arg("sigmaa")))
          .add_property("e_obs_sq", make_getter(&w_t::e_obs_sq, return_value_policy<return_by_value>()))
          .add_property("sig_e_obs_sq", make_getter(&w_t::sig_e_obs_sq, return_value_policy<return_by_value>()))
          .add_property("null_log_z", make_getter(&w_t::null_log_z, return_value_policy<return_by_value>()))
          .add_property("force_exact", make_getter(&w_t::force_exact, return_value_policy<return_by_value>()))
          .add_property("dsqr", make_getter(&w_t::dsqr, return_value_policy<return_by_value>()))
          .def_readonly("rice_kappa", &w_t::rice_kappa)
        ;
      }
      {
        typedef llgi_exact::french_wilson_inverse w_t;
        class_<w_t>("llgi_french_wilson_inverse", no_init)
          .def(init<
            af::const_ref<double> const&,
            af::const_ref<double> const&,
            af::const_ref<double> const&,
            af::const_ref<bool> const&,
            optional<double> >((
              arg("f"),
              arg("sigf"),
              arg("mean_intensity"),
              arg("centric_flags"),
              arg("h_min")=-6.0)))
          .add_property("i_obs", make_getter(&w_t::i_obs, return_value_policy<return_by_value>()))
          .add_property("sig_i_obs", make_getter(&w_t::sig_i_obs, return_value_policy<return_by_value>()))
          .add_property("h", make_getter(&w_t::h, return_value_policy<return_by_value>()))
          .add_property("valid", make_getter(&w_t::valid, return_value_policy<return_by_value>()))
          .add_property("prior_dominated", make_getter(&w_t::prior_dominated, return_value_policy<return_by_value>()))
        ;
      }
      {
        typedef llgi_exact::rice_moments_many w_t;
        class_<w_t>("llgi_rice_moments", no_init)
          .def(init<
            af::const_ref<double> const&,
            af::const_ref<double> const&,
            af::const_ref<bool> const&>((
              arg("e_obs_sq"),
              arg("sig_e_obs_sq"),
              arg("centric_flags"))))
          .add_property("valid", make_getter(&w_t::valid, return_value_policy<return_by_value>()))
          .add_property("mu2", make_getter(&w_t::mu2, return_value_policy<return_by_value>()))
          .add_property("mu4", make_getter(&w_t::mu4, return_value_policy<return_by_value>()))
          .add_property("dsqr", make_getter(&w_t::dsqr, return_value_policy<return_by_value>()))
          .add_property("eeff", make_getter(&w_t::eeff, return_value_policy<return_by_value>()))
        ;
      }
    }
  };

  struct r_factor_wrappers
  {
    //typedef r_factor w_t;

    static void
    wrap()
    {
      using namespace boost::python;

      class_<cctbx::xray::targets::r_factor<> >("r_factor")
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<std::complex<double> > const& >((arg("fo"),arg("fc"))))
      .def("value", &cctbx::xray::targets::r_factor<>::value)
      .def("scale_ls", &cctbx::xray::targets::r_factor<>::scale_ls)
      .def("scale_r", &cctbx::xray::targets::r_factor<>::scale_r)
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
    targets::boost_python::correlation_wrappers::wrap();
    targets::boost_python::least_squares_residual_wrappers<
      cctbx::xray::targets::f_calc_modulus>::wrap(
      "targets_least_squares_residual");
    targets::boost_python::least_squares_residual_wrappers<
      cctbx::xray::targets::f_calc_modulus_square>::wrap(
      "targets_least_squares_residual_for_intensity");
    targets::boost_python::mlf_wrappers::wrap();
    targets::boost_python::mli_wrappers::wrap();
    targets::boost_python::mlhl_wrappers::wrap();
    targets::boost_python::llgi_wrappers::wrap();
    targets::boost_python::llgi_sigmaa_scatfrac_wrappers::wrap();
    targets::boost_python::llgi_e_sigmaa_wrappers::wrap();
    targets::boost_python::llgi_exact_wrappers::wrap();
    targets::boost_python::r_factor_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
