#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/math/mean_and_variance.h>
#include <xfel/metrology/legacy_scale/parameters.h>
#include <xfel/metrology/legacy_scale/bandpass_gaussian.h>
#include <xfel/metrology/legacy_scale/vector_collection.h>
#include <xfel/metrology/legacy_scale/a_g_conversion.h>
#include <xfel/metrology/legacy_scale/quadrants.h>

#include <vector>
#include <map>

using namespace boost::python;
namespace xfel{
namespace boost_python { namespace {

  void
  xfel_scale_init_module() {
    using namespace boost::python;

    typedef return_value_policy<return_by_value> rbv;
    class_<xfel_legacy::parameter::parameter_array>(
      "parameter_array",no_init)
      .add_property("x",
       make_getter(&xfel_legacy::parameter::parameter_array::parameters, rbv()))
      .add_property("gradients",
       make_getter(&xfel_legacy::parameter::parameter_array::gradients,rbv()))
      .add_property("curvatures",
       make_getter(&xfel_legacy::parameter::parameter_array::curvatures,rbv()))
    ;
    class_<xfel_legacy::parameter::organizer_base>(
      "organizer_base",init<>())
      .def("register",&xfel_legacy::parameter::organizer_base::register_array,
        (arg("tag"),arg("ndata"),arg("kdim")=1,arg("data")))
      .def("register_local",&xfel_legacy::parameter::organizer_base::register_local_array,
        (arg("tag"),arg("ndata"),arg("kdim")=1,arg("data")))
      .def("as_x_array",&xfel_legacy::parameter::organizer_base::as_x_array)
      .def("get_gradient_array",
           (xfel_legacy::farray(xfel_legacy::parameter::organizer_base::*)()const)
           &xfel_legacy::parameter::organizer_base::get_gradient_array)
      .def("set_gradient_array",&xfel_legacy::parameter::organizer_base::set_gradient_array)
      .def("get_curvature_array",&xfel_legacy::parameter::organizer_base::get_curvature_array)
      .def("from_x_array",&xfel_legacy::parameter::organizer_base::from_x_array)
      .def("initialize_gradients_curvatures",
         &xfel_legacy::parameter::organizer_base::initialize_gradients_curvatures)
      .def("rezero_gradients_curvatures",
         &xfel_legacy::parameter::organizer_base::rezero_gradients_curvatures)
    ;
    class_<xfel_legacy::algorithm::mark5_iteration,
           bases<xfel_legacy::parameter::organizer_base> >(
      "mark5_iteration",init<>())
      .def("set_refined_origins_to_c",
         &xfel_legacy::algorithm::mark5_iteration::set_refined_origins_to_c)
      .def("compute_target",&xfel_legacy::algorithm::mark5_iteration::compute_target,
        (arg_("tox"),arg_("toy"),
         arg_("spotcx"),arg_("spotcy"),
         arg_("spotfx"),arg_("spotfy"),
         arg_("master_tiles"), arg_("frames"),
         arg_("part_distance")))
      .def("compute_functional_only",
        &xfel_legacy::algorithm::mark5_iteration::compute_functional_only,
        (arg_("tox"),arg_("toy"),
         arg_("spotcx"),arg_("spotcy"),
         arg_("spotfx"),arg_("spotfy"),
         arg_("master_tiles"), arg_("frames"),
         arg_("part_distance")))
      .add_property("model_calcx",
         make_getter(&xfel_legacy::algorithm::mark5_iteration::model_calcx, rbv()))
      .add_property("model_calcy",
         make_getter(&xfel_legacy::algorithm::mark5_iteration::model_calcy, rbv()))
      .def("set_vector_collection",
         &xfel_legacy::algorithm::mark5_iteration::set_vector_collection)
      .def("uncorrected_detector_to_laboratory_frame",
        &xfel_legacy::algorithm::mark5_iteration::uncorrected_detector_to_laboratory_frame,
        (arg_("tox"),arg_("toy"),
         arg_("spotfx"),arg_("spotfy"),
         arg_("master_tiles")
        ))
    ;

    class_<xfel_legacy::parameter::vector_array >(
      "vector_array",init<>())
      .add_property("gradients",
       make_getter(&xfel_legacy::parameter::vector_array::gradients,rbv()))
    ;

    class_<xfel_legacy::parameter::vector_collection >(
      "vector_collection",init<>())
      .def(init<xfel_legacy::iarray,xfel_legacy::marray>())
      .def("collect_vector_information",
        &xfel_legacy::parameter::vector_collection::collect_vector_information)
      .def("register",
        &xfel_legacy::parameter::vector_collection::register_tag,(arg_("tag")))
    ;

    class_<xfel_legacy::parameter::streak_parameters >(
      "streak_parameters",no_init)
      .add_property("rotax",make_getter(
         &xfel_legacy::parameter::streak_parameters::rotax, rbv()))
      .add_property("position",make_getter(
         &xfel_legacy::parameter::streak_parameters::position, rbv()))
      .add_property("rotax_excursion_rad",make_getter(
         &xfel_legacy::parameter::streak_parameters::rotax_excursion_rad, rbv()))
     ;

    class_<xfel_legacy::parameter::bandpass_gaussian>("bandpass_gaussian",
      init<rstbx::bandpass::parameters_bp3 const&>(arg("parameters")))
      .def("set_active_areas", &xfel_legacy::parameter::bandpass_gaussian::set_active_areas)
      .def("set_sensor_model", &xfel_legacy::parameter::bandpass_gaussian::set_sensor_model,(
         arg("thickness_mm"), arg("mu_rho"), arg("signal_penetration")))
      .def("picture_fast_slow_force",
         &xfel_legacy::parameter::bandpass_gaussian::picture_fast_slow_force)
      .def("gaussian_fast_slow",
         &xfel_legacy::parameter::bandpass_gaussian::gaussian_fast_slow)
      .add_property("hi_E_limit",make_getter(
         &xfel_legacy::parameter::bandpass_gaussian::hi_E_limit, rbv()))
      .add_property("lo_E_limit",make_getter(
         &xfel_legacy::parameter::bandpass_gaussian::lo_E_limit, rbv()))
      .add_property("mean_position",make_getter(
         &xfel_legacy::parameter::bandpass_gaussian::mean_position, rbv()))
      .add_property("calc_radial_length",make_getter(
         &xfel_legacy::parameter::bandpass_gaussian::calc_radial_length, rbv()))
      .add_property("part_distance",make_getter(
         &xfel_legacy::parameter::bandpass_gaussian::part_distance, rbv()))
      .add_property("observed_flag",make_getter(
         &xfel_legacy::parameter::bandpass_gaussian::observed_flag, rbv()))
      .def("set_subpixel", &xfel_legacy::parameter::bandpass_gaussian::set_subpixel)
      .def("set_mosaicity", &xfel_legacy::parameter::bandpass_gaussian::set_mosaicity)
      .def("set_domain_size", &xfel_legacy::parameter::bandpass_gaussian::set_domain_size)
      .def("set_bandpass", &xfel_legacy::parameter::bandpass_gaussian::set_bandpass)
      .def("set_orientation", &xfel_legacy::parameter::bandpass_gaussian::set_orientation)
      .def("set_detector_origin",
         &xfel_legacy::parameter::bandpass_gaussian::set_detector_origin)
      .def("set_distance", &xfel_legacy::parameter::bandpass_gaussian::set_distance)
      .def("set_vector_output_pointers",
         &xfel_legacy::parameter::bandpass_gaussian::set_vector_output_pointers,
         (arg_("vector_collection"),arg_("frame_id")))
      .def("measure_bandpass_and_mosaic_parameters", &
      xfel_legacy::parameter::bandpass_gaussian::measure_bandpass_and_mosaic_parameters
         ,(arg("radial"), arg("azimut"), arg("domain_sz_inv_ang"),
           arg("lab_frame_obs")
          ))
      .add_property("wavelength_fit_ang",make_getter(
         &xfel_legacy::parameter::bandpass_gaussian::wavelength_fit_ang, rbv()))
      .add_property("mosaicity_fit_rad",make_getter(
         &xfel_legacy::parameter::bandpass_gaussian::mosaicity_fit_rad, rbv()))
      .add_property("wavelength_fit_ang_sigma",make_getter(
         &xfel_legacy::parameter::bandpass_gaussian::wavelength_fit_ang_sigma, rbv()))
      .add_property("mosaicity_fit_rad_sigma",make_getter(
         &xfel_legacy::parameter::bandpass_gaussian::mosaicity_fit_rad_sigma, rbv()))
      .def("simple_forward_calculation_spot_position",
           &xfel_legacy::parameter::bandpass_gaussian::simple_forward_calculation_spot_position,
           (arg_("wavelength"), arg_("observation_no")))
      .def("simple_part_excursion_part_rotxy",
           &xfel_legacy::parameter::bandpass_gaussian::simple_part_excursion_part_rotxy,
           (arg_("wavelength"), arg_("observation_no"),
            arg_("dA_drotxy")))
    ;

    //class_<cctbx::agconvert::AG>("AGconvert", init<>())
    //  .def("forward", &cctbx::agconvert::AG::forward)
    //  .def("validate_and_setG", &cctbx::agconvert::AG::validate_and_setG)
    //  .def("back_as_orientation", &cctbx::agconvert::AG::back_as_orientation)
    //  .def("back", &cctbx::agconvert::AG::back)
    //  .add_property("G",make_getter(&cctbx::agconvert::AG::G, rbv()))
    //;
    def("best_fit_limit",xfel_legacy::parameter::best_fit_limit);
    def("quadrant_self_correlation",xfel_legacy::qsc);
  }

}
}} // namespace xfel::boost_python::<anonymous>

BOOST_PYTHON_MODULE(xfel_legacy_scale_ext)
{
  xfel::boost_python::xfel_scale_init_module();

}
