#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <cctbx/xray/dispersion_radial.h>

namespace cctbx { namespace xray { namespace boost_python {

namespace {
  template <typename FloatType>
  struct dispersion_radial_correction_wrapper {
    typedef dispersion_radial_correction<FloatType> wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt>("dispersion_radial_correction", no_init)
        .def(init<af::shared<int> const &, int, int, bool, FloatType,
                  FloatType, bool, FloatType>
             ((arg("group_of_scatterer"),
               arg("n_groups"),
               arg("n_terms"),
               arg("grad")=true,
               arg("s_max")=0,
               arg("r_min")=0.1,
               arg("in_cos_two_theta")=true,
               arg("wavelength")=0)))
        .def_readonly("n_terms", &wt::n_terms)
        .def_readonly("n_groups", &wt::n_groups)
        .def_readwrite("s_max", &wt::s_max)
        .def_readwrite("r_min", &wt::r_min)
        .def_readonly("in_cos_two_theta", &wt::in_cos_two_theta)
        .def_readonly("wavelength", &wt::wavelength)
        .def("variable", &wt::variable, (arg("d_star_sq")))
        .def("validate", &wt::validate, (arg("n_samples")=64))
        .def_readwrite("grad", &wt::grad)
        .def_readwrite("grad_index", &wt::grad_index)
        .add_property("group_of_scatterer",
                      make_getter(&wt::group_of_scatterer,
                                  return_value_policy<return_by_value>()))
        .add_property("coefficients",
                      make_getter(&wt::coefficients,
                                  return_value_policy<return_by_value>()))
        .add_property("n_param", &wt::n_param)
        // the array itself, not get_gradients()'s const_ref, which has no
        // to-python converter
        .add_property("gradients",
                      make_getter(&wt::grad_obs,
                                  return_value_policy<return_by_value>()))
        .def("R_at", &wt::R_at, (arg("d_star_sq"), arg("group")))
        .def("fork", &wt::fork)
        ;
      register_ptr_to_python<boost::shared_ptr<wt> >();
    }
  };

} // namespace anonymous

  void wrap_dispersion_radial_correction() {
    dispersion_radial_correction_wrapper<double>::wrap();
  }

}}} // namespace cctbx::xray::boost_python
