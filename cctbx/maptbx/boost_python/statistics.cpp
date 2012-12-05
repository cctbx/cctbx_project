#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/statistics.h>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

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
      typedef mem_iteration<> w_t;
      using namespace boost::python;
      class_<w_t>("mem_iteration", no_init)
        .def(init<af::ref<double, af::c_grid<3> > const&,
                  af::ref<double, af::c_grid<3> > const&,
                  af::ref<double, af::c_grid<3> >,
                  double,
                  af::tiny<int, 3> const&,
                  double,
                  double,
                  bool >())
        .def("tp", &w_t::tp)
        .def("scale", &w_t::scale)
        .def("z", &w_t::z)
        .def("hw", &w_t::hw)
        .def("hn", &w_t::hn)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_statistics()
  {
    statistics_wrappers::wrap();
    more_statistics_wrappers::wrap();
    even_more_statistics_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
