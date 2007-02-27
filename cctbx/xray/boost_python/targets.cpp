#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/targets.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace xray { namespace targets { namespace boost_python {

namespace {

  struct ls_with_scale_wrappers
  {
    typedef ls_with_scale w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("ls_with_scale", no_init)
        .def(init<
          bool,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          af::const_ref< std::complex<double> > const&,
          bool const&,
          double>((
            arg_("apply_scale_to_f_calc"),
            arg_("f_obs"),
            arg_("weights"),
            arg_("f_calc"),
            arg_("compute_derivatives"),
            arg_("scale_factor"))))
        .def("apply_scale_to_f_calc", &w_t::apply_scale_to_f_calc)
        .def("target", &w_t::target)
        .def("derivatives", &w_t::derivatives, ccr())
        .def("scale_factor", &w_t::scale_factor)
      ;
    }
  };

  struct least_squares_residual_wrappers
  {
    typedef least_squares_residual<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("targets_least_squares_residual",
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
    typedef maximum_likelihood_criterion<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("targets_maximum_likelihood_criterion",
             init<af::const_ref<double> const&,
                  af::const_ref<std::complex<double> > const&,
                  af::const_ref<double> const&,
                  af::const_ref<double> const&,
                  double const&,
                  af::const_ref<int> const&,
                  af::const_ref<bool> const&,
                  bool const& >())
        .def("target", &w_t::target)
        .def("derivatives", &w_t::derivatives)
      ;
    }
  };

  struct maximum_likelihood_criterion_hl_wrappers
  {
    typedef maximum_likelihood_criterion_hl<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("targets_maximum_likelihood_criterion_hl",
             init<af::const_ref<double> const&,
                  af::const_ref<std::complex<double> > const&,
                  af::const_ref<double> const&,
                  af::const_ref<double> const&,
                  af::const_ref<int> const&,
                  af::const_ref<bool> const&,
                  bool const&,
                  af::const_ref<cctbx::hendrickson_lattman<double> > const&,
                  double const& >())
        .def("target", &w_t::target)
        .def("targets", &w_t::targets)
        .def("derivatives", &w_t::derivatives)
      ;
    }
  };

} // namespace <anoymous>

}} // namespace targets::boost_python

namespace boost_python {

  void wrap_targets()
  {
    targets::boost_python::ls_with_scale_wrappers::wrap();
    targets::boost_python::least_squares_residual_wrappers::wrap();
    targets::boost_python::intensity_correlation_wrappers::wrap();
    targets::boost_python::maximum_likelihood_criterion_wrappers::wrap();
    targets::boost_python::maximum_likelihood_criterion_hl_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
