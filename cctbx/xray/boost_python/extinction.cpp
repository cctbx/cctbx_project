#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <cctbx/xray/extinction.h>

namespace cctbx { namespace xray { namespace boost_python {

namespace {
  template <typename FloatType>
  struct extinction_correction_wrapper {
    typedef extinction_correction<FloatType> wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt, boost::noncopyable>("extinction_correction", no_init)
        .def("compute", &wt::compute,
          (arg("index"),
           arg("fc"),
           arg("gradient"),
           arg("compute_gradient")));
    }
  };

  template <typename FloatType>
  struct dummy_extinction_correction_wrapper {
    typedef dummy_extinction_correction<FloatType> wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
        bases<extinction_correction<FloatType> > >
          ("dummy_extinction_correction", no_init)
        .def(init<>())
        .add_property("grad", &wt::grad_value)
        ;
    }
  };

  template <typename FloatType>
  struct shelx_extinction_correction_wrapper {
    typedef shelx_extinction_correction<FloatType> wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
        bases<extinction_correction<FloatType> > >
          ("shelx_extinction_correction", no_init)
        .def(init<uctbx::unit_cell const &,
                  FloatType,
                  FloatType>
             ((arg("unit_cell"),
               arg("wavelength"),
               arg("value"))))
        .def_readwrite("value", &wt::value)
        .def_readwrite("grad_index", &wt::grad_index)
        .def_readwrite("grad", &wt::grad)
        ;
    }
  };

} // namespace anonymous

  void wrap_extinction_correction() {
    extinction_correction_wrapper<double>::wrap();
    dummy_extinction_correction_wrapper<double>::wrap();
    shelx_extinction_correction_wrapper<double>::wrap();
  }

}}} // namespace cctbx::xray::boost_python
