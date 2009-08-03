#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/math/basic_statistics.h>

#include <boost/python/class.hpp>

namespace scitbx { namespace af { namespace boost_python {

  template <typename FloatType>
  struct median_functor_wrapper
  {
    typedef math::median_functor wt;

    static
    FloatType
    call(
      wt& O,
      af::ref<FloatType> const& data) { return O(data); }

    static
    math::median_statistics<FloatType>
    dispersion(
      wt& O,
      af::ref<FloatType> const &data) { return O.dispersion(data); }

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<>())
        .def(init<wt::random_number_engine_t::result_type>(arg("seed")))
/* Intel C++ 9.1 does not support this syntax:
        .def("__call__", &wt::operator()<FloatType>, arg("data"))
        .def("dispersion", &wt::dispersion<FloatType>, arg("data"))
   Replaced with thin wrappers, for simplicity on all platforms.
 */
        .def("__call__", call, arg("data"))
        .def("dispersion", dispersion, arg("data"))
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
