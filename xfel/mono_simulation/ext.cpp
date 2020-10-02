#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/math/mean_and_variance.h>
#include <xfel/mono_simulation/parameters.h>
#include <xfel/mono_simulation/bandpass_gaussian.h>
#include <xfel/mono_simulation/vector_collection.h>


#include <vector>
#include <map>

using namespace boost::python;
namespace xfel{
namespace boost_python { namespace {

  void
  xfel_mono_sim_init_module() {
    using namespace boost::python;

    typedef return_value_policy<return_by_value> rbv;
    class_<xfel::parameter::parameter_array>(
      "parameter_array",no_init)
      .add_property("x",
       make_getter(&xfel::parameter::parameter_array::parameters, rbv()))
      .add_property("gradients",
       make_getter(&xfel::parameter::parameter_array::gradients,rbv()))
      .add_property("curvatures",
       make_getter(&xfel::parameter::parameter_array::curvatures,rbv()))
    ;
    class_<xfel::parameter::organizer_base>(
      "organizer_base",init<>())
      .def("register",&xfel::parameter::organizer_base::register_array,
        (arg("tag"),arg("ndata"),arg("kdim")=1,arg("data")))
      .def("register_local",&xfel::parameter::organizer_base::register_local_array,
        (arg("tag"),arg("ndata"),arg("kdim")=1,arg("data")))
      .def("as_x_array",&xfel::parameter::organizer_base::as_x_array)
      .def("get_gradient_array",
           (xfel::farray(xfel::parameter::organizer_base::*)()const)
           &xfel::parameter::organizer_base::get_gradient_array)
      .def("set_gradient_array",&xfel::parameter::organizer_base::set_gradient_array)
      .def("get_curvature_array",&xfel::parameter::organizer_base::get_curvature_array)
      .def("from_x_array",&xfel::parameter::organizer_base::from_x_array)
      .def("initialize_gradients_curvatures",
         &xfel::parameter::organizer_base::initialize_gradients_curvatures)
      .def("rezero_gradients_curvatures",
         &xfel::parameter::organizer_base::rezero_gradients_curvatures)
    ;

    class_<xfel::parameter::vector_array >(
      "vector_array",init<>())
      .add_property("gradients",
       make_getter(&xfel::parameter::vector_array::gradients,rbv()))
    ;

    class_<xfel::parameter::vector_collection >(
      "vector_collection",init<>())
      .def(init<xfel::iarray,xfel::marray>())
      .def("collect_vector_information",
        &xfel::parameter::vector_collection::collect_vector_information)
      .def("register",
        &xfel::parameter::vector_collection::register_tag,(arg_("tag")))
    ;

    class_<xfel::parameter::streak_parameters >(
      "streak_parameters",no_init)
      .add_property("rotax",make_getter(
         &xfel::parameter::streak_parameters::rotax, rbv()))
      .add_property("position",make_getter(
         &xfel::parameter::streak_parameters::position, rbv()))
      .add_property("position_to_fictitious",make_getter(
         &xfel::parameter::streak_parameters::position_to_fictitious, rbv()))
      .add_property("rotax_excursion_rad",make_getter(
         &xfel::parameter::streak_parameters::rotax_excursion_rad, rbv()))
      .add_property("rotax_excursion_rad_pvr",make_getter(
         &xfel::parameter::streak_parameters::rotax_excursion_rad_pvr, rbv()))
     ;

    class_<xfel::parameter::bandpass_gaussian>("bandpass_gaussian",
      init<rstbx::bandpass::parameters_bp3 const&>(arg("parameters")))
      .def("set_active_areas", &xfel::parameter::bandpass_gaussian::set_active_areas)
      .def("set_sensor_model", &xfel::parameter::bandpass_gaussian::set_sensor_model,(
         arg("thickness_mm"), arg("mu_rho"), arg("signal_penetration")))
      .def("gaussian_fast_slow",
         &xfel::parameter::bandpass_gaussian::gaussian_fast_slow)
      .add_property("hi_E_limit",make_getter(
         &xfel::parameter::bandpass_gaussian::hi_E_limit, rbv()))
      .add_property("lo_E_limit",make_getter(
         &xfel::parameter::bandpass_gaussian::lo_E_limit, rbv()))
      .add_property("mean_position",make_getter(
         &xfel::parameter::bandpass_gaussian::mean_position, rbv()))
      .add_property("calc_radial_length",make_getter(
         &xfel::parameter::bandpass_gaussian::calc_radial_length, rbv()))
      .add_property("part_distance",make_getter(
         &xfel::parameter::bandpass_gaussian::part_distance, rbv()))
      .add_property("observed_flag",make_getter(
         &xfel::parameter::bandpass_gaussian::observed_flag, rbv()))
      .def("set_subpixel", &xfel::parameter::bandpass_gaussian::set_subpixel,(
           arg("translations"), arg("rotations_deg")))
      .def("set_mosaicity", &xfel::parameter::bandpass_gaussian::set_mosaicity)
      .def("set_domain_size", &xfel::parameter::bandpass_gaussian::set_domain_size)
      .def("set_bandpass", &xfel::parameter::bandpass_gaussian::set_bandpass)
      .def("set_orientation", &xfel::parameter::bandpass_gaussian::set_orientation)
      .def("set_detector_origin",
         &xfel::parameter::bandpass_gaussian::set_detector_origin)
      .def("set_distance", &xfel::parameter::bandpass_gaussian::set_distance)
      .def("set_vector_output_pointers",
         &xfel::parameter::bandpass_gaussian::set_vector_output_pointers,
         (arg_("vector_collection"),arg_("frame_id")))
      .def("measure_bandpass_and_mosaic_parameters", &
      xfel::parameter::bandpass_gaussian::measure_bandpass_and_mosaic_parameters
         ,(arg("radial"), arg("azimut"), arg("domain_sz_inv_ang"),
           arg("lab_frame_obs")
          ))
      .add_property("wavelength_fit_ang",make_getter(
         &xfel::parameter::bandpass_gaussian::wavelength_fit_ang, rbv()))
      .add_property("mosaicity_fit_rad",make_getter(
         &xfel::parameter::bandpass_gaussian::mosaicity_fit_rad, rbv()))
      .add_property("wavelength_fit_ang_sigma",make_getter(
         &xfel::parameter::bandpass_gaussian::wavelength_fit_ang_sigma, rbv()))
      .add_property("mosaicity_fit_rad_sigma",make_getter(
         &xfel::parameter::bandpass_gaussian::mosaicity_fit_rad_sigma, rbv()))
      .def("simple_forward_calculation_spot_position",
           &xfel::parameter::bandpass_gaussian::simple_forward_calculation_spot_position,
           (arg_("wavelength"), arg_("observation_no")))
      .def("simple_part_excursion_part_rotxy",
           &xfel::parameter::bandpass_gaussian::simple_part_excursion_part_rotxy,
           (arg_("wavelength"), arg_("observation_no"),
            arg_("dA_drotxy")))
    ;

    def("best_fit_limit",xfel::parameter::best_fit_limit);
  }

}
}} // namespace xfel::boost_python::<anonymous>

BOOST_PYTHON_MODULE(xfel_mono_sim_ext)
{
  xfel::boost_python::xfel_mono_sim_init_module();

}
