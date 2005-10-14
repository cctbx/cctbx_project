#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/mean_and_variance.h>
#include <boost/python/class.hpp>

namespace scitbx { namespace af { namespace boost_python { namespace {

  struct mean_and_variance_wrappers
  {
    typedef math::mean_and_variance<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("mean_and_variance", no_init)
        .def(init<af::const_ref<double> const&>())
        .def(init<af::const_ref<double> const&,
                  af::const_ref<double> const&>())
        .def("have_weights", &w_t::have_weights)
        .def("mean", &w_t::mean)
        .def("gsl_stats_wvariance", &w_t::gsl_stats_wvariance)
        .def("gsl_stats_wsd", &w_t::gsl_stats_wsd)
        .def("standard_error_of_mean_calculated_from_sample_weights",
          &w_t::standard_error_of_mean_calculated_from_sample_weights)
        .def("unweighted_sample_variance", &w_t::unweighted_sample_variance)
        .def("unweighted_sample_standard_deviation",
          &w_t::unweighted_sample_standard_deviation)
        .def("unweighted_standard_error_of_mean",
          &w_t::unweighted_standard_error_of_mean)
        .def("sum_weights", &w_t::sum_weights)
        .def("sum_weights_sq", &w_t::sum_weights_sq)
        .def("sum_weights_values", &w_t::sum_weights_values)
        .def("sum_weights_delta_sq", &w_t::sum_weights_delta_sq)
      ;
    }
  };

} // namespace <anonymous>

  void wrap_flex_mean_and_variance()
  {
    mean_and_variance_wrappers::wrap();
  }

}}} // namespace scitbx::af::boost_python
