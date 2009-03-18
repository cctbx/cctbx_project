#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/module.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <cctbx/cumulative_intensity_distribution.h>

namespace cctbx { namespace boost_python { namespace statistics_ext {

  template <typename FloatType>
  struct cumulative_intensity_wrapper
  {
    typedef cumulative_intensity<FloatType> wt;

    static void wrap() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      typedef return_internal_reference<> rir;

      class_<wt>("cumulative_intensity_core", no_init)
        .def(init<af::const_ref<FloatType> const &,
                  af::const_ref<FloatType> const &,
                  af::const_ref<FloatType> const &,
                  af::const_ref<FloatType> const &,
                  af::shared<miller::index<> > const &>(
              (arg("data"),
               arg("d_spacings"),
               arg("mean_data"),
               arg("bin_d_max"),
               arg("indices"))))
        .def("x", make_getter(&wt::x_, rbv))
        .def("y", make_getter(&wt::y_, rbv))
        ;
    }
  };

  void wrap_cumulative_intensity() {
    cumulative_intensity_wrapper<double>::wrap();
  }

  void init_module()
  {
    using namespace boost::python;

    wrap_cumulative_intensity();
  }
}}} // cctbx::boostpython::statistics_ext

BOOST_PYTHON_MODULE(cctbx_statistics_ext)
{
  cctbx::boost_python::statistics_ext::init_module();
}
