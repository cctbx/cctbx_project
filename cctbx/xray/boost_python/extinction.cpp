#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <cctbx/xray/extinction.h>

namespace cctbx { namespace xray { namespace boost_python {

namespace {
  template <typename FloatType>
  struct fc_correction_wrapper {
    typedef fc_correction<FloatType> wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt, boost::noncopyable>("fc_correction", no_init)
        .def("compute", &wt::compute, (
          arg("index"),
          arg("fc_sq"),
          arg("compute_gradient")))
        .add_property("grad_fc_multiplier", &wt::get_grad_Fc_multiplier)
        .add_property("grad_index", &wt::get_grad_index)
        .add_property("gradients", &wt::get_gradients)
        .def_readwrite("grad", &wt::grad)
        ;
    }
  };


  template <typename FloatType>
  struct dummy_fc_correction_wrapper {
    typedef dummy_fc_correction<FloatType> wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
        bases<fc_correction<FloatType> > >
          ("dummy_fc_correction", no_init)
        .def(init<>())
        ;
    }
  };

  template <typename FloatType>
  struct shelx_extinction_correction_wrapper {
    typedef shelx_extinction_correction<FloatType> wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
        bases<fc_correction<FloatType> > >
          ("shelx_extinction_correction", no_init)
        .def(init<uctbx::unit_cell const &,
                  FloatType,
                  FloatType>
             ((arg("unit_cell"),
               arg("wavelength"),
               arg("value"))))
        .def_readwrite("value", &wt::value)
        .def_readwrite("grad_index", &wt::grad_index)
        ;
    }
  };

  template <typename FloatType>
  struct shelx_SWAT_correction_wrapper {
    typedef shelx_SWAT_correction<FloatType> wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
        bases<fc_correction<FloatType> > >
        ("shelx_SWAT_correction", no_init)
        .def(init<uctbx::unit_cell const&,
          FloatType,
          FloatType>
          ((arg("unit_cell"),
            arg("g"),
            arg("U"))))
        .add_property("g", &wt::get_g, &wt::set_g)
        .add_property("U", &wt::get_U, &wt::set_U)
        .def_readwrite("grad_index", &wt::grad_index)
        ;
    }
  };

} // namespace anonymous

  void wrap_extinction_correction() {
    fc_correction_wrapper<double>::wrap();
    dummy_fc_correction_wrapper<double>::wrap();
    shelx_extinction_correction_wrapper<double>::wrap();
    shelx_SWAT_correction_wrapper<double>::wrap();
  }

}}} // namespace cctbx::xray::boost_python
