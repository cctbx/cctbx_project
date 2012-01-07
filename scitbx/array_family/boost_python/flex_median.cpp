#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/math/basic_statistics.h>

#include <boost/python/class.hpp>
#include <vector>

namespace scitbx { namespace af { namespace boost_python {

  template <typename FloatType>
  struct median_functor_wrapper
  {
    typedef math::median_functor wt;

    static
    FloatType
    call(
      wt& O,
      af::const_ref<FloatType> const& data)
    {
      std::vector<FloatType> buffer(data.begin(), data.end());
      return O(af::make_ref(buffer));
    }

    static
    math::median_statistics<FloatType>
    dispersion(
      wt& O,
      af::const_ref<FloatType> const &data)
    {
      std::vector<FloatType> buffer(data.begin(), data.end());
      return O.dispersion(af::make_ref(buffer));
    }

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<>())
        .def(init<wt::random_number_engine_t::result_type>(arg("seed")))
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
