#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/math/basic_statistics.h>

#include <boost/python/class.hpp>

namespace scitbx { namespace af { namespace boost_python {

  template <typename FloatType>
  struct median_functor_wrapper
  {
    typedef math::median_functor wt;

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<>())
        .def(init<wt::random_number_engine_t::result_type>(arg("seed")))
        .def("__call__", &wt::operator()<FloatType>, arg("data"))
        .def("dispersion", &wt::dispersion<FloatType>, arg("data"))
      ;
    }
  };


  template <typename FloatType>
  struct median_statistics_wrapper
  {
    typedef math::median_statistics<FloatType> wt;

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def_readonly("median", &wt::median)
        .def_readonly("median_absolute_deviation",
                      &wt::median_absolute_deviation)
      ;
    }
  };

  void wrap_flex_median_statistics() {
    median_statistics_wrapper<double>::wrap("median_statistics");
    median_functor_wrapper<double>::wrap("median_functor");
  }

}}}
