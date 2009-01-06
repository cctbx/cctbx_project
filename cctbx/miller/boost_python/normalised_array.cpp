#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <scitbx/array_family/boost_python/shared_wrapper.h>


#include <cctbx/miller/normalised_array.h>

namespace cctbx { namespace miller { namespace boost_python {

  template <typename FloatType>
  struct normalised_array_wrapper
  {
    typedef normalised_array<FloatType> wt;

    static void wrap() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      typedef return_internal_reference<> rir;

      af::boost_python::shared_wrapper<
        typename wt::form_factor_t, rir>::wrap(
        "shared_gaussian_form_factors");

      class_<wt>("normalised_array", no_init)
        .def(init<af::const_ref<typename wt::form_factor_t> const &,
                  af::const_ref<FloatType> const &,
                  FloatType,
                  FloatType,
                  uctbx::unit_cell const &,
                  sgtbx::space_group const &,
                  af::shared<index<> > const &,
                  af::shared<FloatType> const &>(
              (arg("form_factors"),
               arg("multiplicities"),
               arg("wilson_intensity_scale_factor"),
               arg("wilson_b"),
               arg("unit_cell"),
               arg("space_group"),
               arg("indices"),
               arg("data"))))
        .def("data", make_getter(&wt::data_out, rbv))
        .def("sum_e_sq_minus_1", make_getter(&wt::sum_e_sq_minus_1))
        .def("n_e_greater_than_2", make_getter(&wt::n_e_greater_than_2))
        ;
    }
  };

  void wrap_normalised_array() {
    normalised_array_wrapper<double>::wrap();
  }
}}} // cctbx::miller::boostpython
