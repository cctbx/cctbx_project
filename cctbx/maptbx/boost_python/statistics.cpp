#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/statistics.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  struct statistics_wrappers
  {
    typedef statistics<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("statistics", no_init)
        .def(init<af::const_ref<float, af::flex_grid<> > const&>())
        .def(init<af::const_ref<double, af::flex_grid<> > const&>())
        .def("min", &w_t::min)
        .def("max", &w_t::max)
        .def("mean", &w_t::mean)
        .def("mean_sq", &w_t::mean_sq)
        .def("sigma", &w_t::sigma)
      ;
    }
  };

  struct more_statistics_wrappers
  {
    typedef more_statistics<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<statistics<> > >("more_statistics", no_init)
        .def(init<af::const_ref<float, af::flex_grid<> > const&>())
        .def(init<af::const_ref<double, af::flex_grid<> > const&>())
        .def("skewness", &w_t::skewness)
        .def("kurtosis", &w_t::kurtosis)
      ;
    }
  };

  struct even_more_statistics_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<update_prior>("update_prior", no_init)
        .def(init<af::const_ref<std::complex<double>, af::flex_grid<> > const&,
                  af::const_ref<std::complex<double>, af::flex_grid<> > const&,
                  af::versa<std::complex<double>, af::flex_grid<> > >((
          arg("fobs"),
          arg("sigf"),
          arg("priorA"))))
        .def_readonly("chi2", &update_prior::chi2)
        .def_readonly("sum", &update_prior::sum);
    }
  };

} // namespace <anoymous>

  void wrap_statistics()
  {
    statistics_wrappers::wrap();
    more_statistics_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
