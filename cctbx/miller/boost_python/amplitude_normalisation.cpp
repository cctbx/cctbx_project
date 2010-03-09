#include <boost/python/class.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <scitbx/array_family/boost_python/shared_wrapper.h>


#include <cctbx/miller/amplitude_normalisation.h>

namespace cctbx { namespace miller { namespace boost_python {

  template <typename FloatType>
  struct amplitude_normalisation_wrapper
  {
    typedef amplitude_normalisation<FloatType> wt;
    typedef typename wt::float_type float_type;

    static void wrap() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      typedef return_internal_reference<> rir;

      af::boost_python::shared_wrapper<
        typename wt::form_factor_t, rir>::wrap(
        "shared_gaussian_form_factors");

      class_<wt>("amplitude_normalisation", no_init)
        .def(init<af::const_ref<typename wt::form_factor_t> const &,
                  af::const_ref<float_type> const &,
                  float_type,
                  float_type,
                  uctbx::unit_cell const &,
                  sgtbx::space_group const &,
                  af::const_ref<index<> > const &>(
              (arg("form_factors"),
               arg("multiplicities"),
               arg("wilson_intensity_scale_factor"),
               arg("wilson_b"),
               arg("unit_cell"),
               arg("space_group"),
               arg("indices"))))
        .add_property("normalisations", make_getter(&wt::normalisations, rbv))
        ;
    }
  };

  void wrap_amplitude_normalisation() {
    amplitude_normalisation_wrapper<double>::wrap();
  }
}}} // cctbx::miller::boostpython
