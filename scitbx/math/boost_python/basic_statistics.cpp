#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/basic_statistics.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

namespace scitbx { namespace math { namespace {

  struct basic_statistics_wrappers
  {
    typedef basic_statistics<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("basic_statistics", no_init)
        .def(init<af::const_ref<double> const&>((arg_("values"))))
        .def_readonly("n", &w_t::n)
        .def_readonly("min", &w_t::min)
        .def_readonly("max", &w_t::max)
        .def_readonly("max_absolute", &w_t::max_absolute)
        .def_readonly("sum", &w_t::sum)
        .def_readonly("mean", &w_t::mean)
        .def_readonly("mean_absolute_deviation_from_mean",
          &w_t::mean_absolute_deviation_from_mean)
        .def_readonly("biased_variance", &w_t::biased_variance)
        .def_readonly("biased_standard_deviation",
          &w_t::biased_standard_deviation)
        .def_readonly("bias_corrected_variance", &w_t::bias_corrected_variance)
        .def_readonly("bias_corrected_standard_deviation",
          &w_t::bias_corrected_standard_deviation)
        .def_readonly("skew", &w_t::skew)
        .def_readonly("kurtosis", &w_t::kurtosis)
        .def_readonly("kurtosis_excess", &w_t::kurtosis_excess)
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_basic_statistics()
  {
    basic_statistics_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
