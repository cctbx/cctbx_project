#include <boost/python/class.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <scitbx/array_family/boost_python/shared_wrapper.h>


#include <cctbx/miller/amplitude_rescaling.h>

namespace cctbx { namespace miller { namespace boost_python {

  template <typename FloatType>
  struct amplitude_rescaling_wrapper
  {
    typedef amplitude_rescaling<FloatType> wt;

    static void wrap() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      typedef return_internal_reference<> rir;

      af::boost_python::shared_wrapper<
        typename wt::form_factor_t, rir>::wrap(
        "shared_gaussian_form_factors");

      enum_<typename wt::rescaling_kind>("amplitude_rescaling_kind")
        .value("normalised", wt::normalised)
        .value("denormalised", wt::denormalised)
        ;
      class_<wt>("amplitude_rescaling", no_init)
        .def(init<int,
                  af::const_ref<typename wt::form_factor_t> const &,
                  af::const_ref<FloatType> const &,
                  FloatType,
                  FloatType,
                  uctbx::unit_cell const &,
                  sgtbx::space_group const &,
                  af::shared<index<> > const &,
                  af::shared<FloatType> const &>(
              (arg("kind"),
               arg("form_factors"),
               arg("multiplicities"),
               arg("wilson_intensity_scale_factor"),
               arg("wilson_b"),
               arg("unit_cell"),
               arg("space_group"),
               arg("indices"),
               arg("data"))))
        .def("data", make_getter(&wt::result, rbv))
        .def("sum_e_sq_minus_1", make_getter(&wt::sum_e_sq_minus_1))
        .def("n_e_greater_than_2", make_getter(&wt::n_e_greater_than_2))
        ;
    }
  };

  void wrap_amplitude_rescaling() {
    amplitude_rescaling_wrapper<double>::wrap();
  }
}}} // cctbx::miller::boostpython
