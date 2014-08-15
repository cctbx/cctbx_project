#include <scitbx/math/weighted_covariance.h>

#include <boost/python/class.hpp>
#include <boost/python/return_arg.hpp>

namespace scitbx { namespace math { namespace boost_python {

  template <typename FloatType>
  struct weighted_covariance_wrapper
  {
    typedef weighted_covariance<FloatType> wt;
    typedef af::const_ref<FloatType> const &const_ref_t;

    #define ADD_PROPERTY(name) add_property(#name, &wt::name)

    static void wrap(char const *name) {
      using namespace boost::python;
      return_self<> rs;
      class_<wt>(name, no_init)
        .def(init<>())
        .def(init<const_ref_t, const_ref_t, const_ref_t>
             ((arg("x"), arg("y"), arg("weights"))))
        .def("accumulate", &wt::accumulate,
             (arg("x"), arg("y"), arg("weight")), rs)
        .ADD_PROPERTY(mean_x)
        .ADD_PROPERTY(mean_y)
        .ADD_PROPERTY(variance_x)
        .ADD_PROPERTY(variance_y)
        .ADD_PROPERTY(covariance_xy)
        .ADD_PROPERTY(correlation)
        ;
    }
  };


  template <typename FloatType>
  struct multivariate_moments_wrapper
  {
    typedef multivariate_moments<FloatType> wt;
    typedef af::const_ref<FloatType> const &const_ref_t;

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name,no_init)
        .def(init<>())
        .def(init<const_ref_t>
              ((arg("weights"))))
        .def("update", &wt::update,
               (arg("data")) )
        .def("mean", &wt::mean)
        .def("variance", &wt::variance)
        .def("vcv_upper_triangle_packed", &wt::vcv_upper_triangle_packed)
        .def("vcv_raw_upper_triangle_packed", &wt::vcv_raw_upper_triangle_packed)
        ;
    }


  };


  void wrap_weighted_covariance() {
    weighted_covariance_wrapper<double>::wrap("weighted_covariance");
    multivariate_moments_wrapper<double>::wrap("multivariate_moments");
  }


}}}
