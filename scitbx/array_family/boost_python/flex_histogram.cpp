#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/histogram.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(
  scitbx::histogram<>)

namespace scitbx { namespace af { namespace boost_python { namespace {

  struct histogram_wrappers
  {
    typedef histogram<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("histogram", no_init)
        .def(init<
          af::const_ref<double> const&,
          std::size_t>((
            arg("data"),
            arg("n_slots")=1000)))
        .def(init<
          af::const_ref<double> const&,
          double const&,
          double const&,
          std::size_t,
          double const&>((
            arg("data"),
            arg("data_min"),
            arg("data_max"),
            arg("n_slots")=1000,
            arg("relative_tolerance")=1e-4)))
        .def(init<
          w_t const&,
          af::const_ref<double> const&,
          double const&>((
            arg("other"),
            arg("data"),
            arg("relative_tolerance")=1e-4)))
        .enable_pickling()
        .def(init<
          double const&,
          double const&,
          double const&,
          af::shared<long> const&,
          std::size_t>()) // intentionally no arg
        .def("data_min", &w_t::data_min)
        .def("data_max", &w_t::data_max)
        .def("slot_width", &w_t::slot_width)
        .def("slot_centers", &w_t::slot_centers)
        .def("slots", &w_t::slots)
        .def("n_out_of_slot_range", &w_t::n_out_of_slot_range)
        .def("get_cutoff", &w_t::get_cutoff, (
          arg("max_points"),
          arg("relative_tolerance")=1e-4))
      ;
    }
  };

} // namespace <anonymous>

  void wrap_flex_histogram()
  {
    histogram_wrappers::wrap();
  }

}}} // namespace scitbx::af::boost_python
