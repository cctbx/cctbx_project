#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python.hpp>
#include <simtbx/diffBragg/src/diffBragg.h>
//#include <boost/python/tuple.hpp>

using namespace boost::python;
namespace simtbx{
namespace diffBragg{
namespace boost_python { namespace {

  static void  set_Nabc(simtbx::nanoBragg::diffBragg& diffBragg, double const& value) {
      //diffBragg.Na = value;
      //diffBragg.Nb = value;
      //diffBragg.Nc = value;
      ///* clear things that might override it */
      //diffBragg.xtal_size_x = -1;
      //diffBragg.xtal_size_y = -1;
      //diffBragg.xtal_size_z = -1;
      //diffBragg.update_oversample();
      diffBragg.set_value(9, value); // FIXME 9 means Ncells parameter, fix this API
  }

  static double get_Nabc(simtbx::nanoBragg::diffBragg const& diffBragg) {
      return diffBragg.Na;
  }

  // change of basis operator see dxtbx.model.crystal.h
  static nanoBragg::mat3 get_O(simtbx::nanoBragg::diffBragg& diffBragg){
    return diffBragg.Omatrix;
  }

  void set_O(simtbx::nanoBragg::diffBragg& diffBragg, nanoBragg::mat3 const& value){
    diffBragg.Omatrix = value;
  }

  static nanoBragg::mat3 get_U(simtbx::nanoBragg::diffBragg& diffBragg){
    return diffBragg.Umatrix;
  }

  void set_U(simtbx::nanoBragg::diffBragg& diffBragg, nanoBragg::mat3 const& value){
    diffBragg.Umatrix = value;
  }

  static nanoBragg::mat3 get_B(simtbx::nanoBragg::diffBragg& diffBragg){
    return diffBragg.Bmatrix;
  }

  void set_B(simtbx::nanoBragg::diffBragg& diffBragg, nanoBragg::mat3 const& value){
      diffBragg.Bmatrix = value;
      diffBragg.a_star[1] = value(0,0);  // initialize Angstrom version, update at the end
      diffBragg.a_star[2] = value(1,0);
      diffBragg.a_star[3] = value(2,0);
      diffBragg.b_star[1] = value(0,1);
      diffBragg.b_star[2] = value(1,1);
      diffBragg.b_star[3] = value(2,1);
      diffBragg.c_star[1] = value(0,2);
      diffBragg.c_star[2] = value(1,2);
      diffBragg.c_star[3] = value(2,2);

      diffBragg.user_cell = 0;
      diffBragg.user_matrix = 1;
      diffBragg.init_cell();
  }

  void set_roi(simtbx::nanoBragg::diffBragg& diffBragg, boost::python::tuple const& tupleoftuples){
      boost::python::tuple frange;
      boost::python::tuple srange;
      frange = extract<boost::python::tuple>(tupleoftuples[0]);
      srange = extract<boost::python::tuple>(tupleoftuples[1]);
      diffBragg.roi_xmin=extract<int>(frange[0]);
      diffBragg.roi_xmax=extract<int>(frange[1]);
      diffBragg.roi_ymin=extract<int>(srange[0]);
      diffBragg.roi_ymax=extract<int>(srange[1]);
       /* fix this up so its not out-of-range */
      diffBragg.init_beamcenter();
      diffBragg.init_raw_pixels_roi();
      diffBragg.initialize_managers();
  }

  /* region of interest in pixels */
  static boost::python::tuple get_roi(simtbx::nanoBragg::diffBragg const& diffBragg) {
      boost::python::tuple frange;
      boost::python::tuple srange;
      frange = boost::python::make_tuple(diffBragg.roi_xmin,diffBragg.roi_xmax);
      srange = boost::python::make_tuple(diffBragg.roi_ymin,diffBragg.roi_ymax);
      return boost::python::make_tuple(frange,srange);
  }

  void diffBragg_init_module() {
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;
    typedef return_internal_reference<> rir;
    class_<simtbx::nanoBragg::diffBragg, bases<simtbx::nanoBragg::nanoBragg> >
            ("diffBragg", no_init)
      /* constructor that takes a dxtbx detector and beam model */
      .def(init<const dxtbx::model::Detector&,
                const dxtbx::model::Beam&,
                int, int>(
        (arg_("detector"),
         arg_("beam"),
         arg_("verbose")=0,
         arg_("panel_id")=0),
        "nanoBragg simulation initialized from dxtbx detector and beam objects"))

      /* function for doing some differentiation */
      .def("add_diffBragg_spots",&simtbx::nanoBragg::diffBragg::add_diffBragg_spots,
        "gives derivitive of average photon count w.r.t. parameter of choice")

      /* Run this before running add_diffBragg_spots */
      .def("vectorize_umats",&simtbx::nanoBragg::diffBragg::vectorize_umats,
        "caches the UMATS")

      .def("initialize_managers",&simtbx::nanoBragg::diffBragg::initialize_managers,
        "sets up the managers")

      .def("refine", &simtbx::nanoBragg::diffBragg::refine,
        "sets the refinement flag")

      .def("get_value", &simtbx::nanoBragg::diffBragg::get_value, "get value of the refinement parameter")

      .def("set_value", &simtbx::nanoBragg::diffBragg::set_value, "set value of the refinement parameter")

      .def("set_ucell_derivative_matrix",  &simtbx::nanoBragg::diffBragg::set_ucell_derivative_matrix, "Boo-ya")

      .def("set_ucell_second_derivative_matrix",  &simtbx::nanoBragg::diffBragg::set_ucell_second_derivative_matrix,
            "Two-ya")

      //.def("get_ucell_derivative_matrix",  &simtbx::nanoBragg::diffBragg::get_ucell_derivative_matrix, "scooby snacks")

      .def("get_derivative_pixels", &simtbx::nanoBragg::diffBragg::get_derivative_pixels,
            "gets the manager raw image containing first derivatives")

      .def("get_second_derivative_pixels", &simtbx::nanoBragg::diffBragg::get_second_derivative_pixels,
            "gets the manager raw image containing second derivatives")

      .def("zero_raw_pixel_rois", &simtbx::nanoBragg::diffBragg::zero_raw_pixel_rois,
           "reallocate the raw image ROIs (they are usually tiny)")

      .def("update_dxtbx_geoms", &simtbx::nanoBragg::diffBragg::update_dxtbx_geoms,
           "update the geometries with new dxtbx models, number of pixels should remain constant")

      .add_property("region_of_interest",
             make_function(&get_roi,rbv()),
             make_function(&set_roi,dcp()),
             "region of interest on detector: fast_min fast_max slow_min slow_max")

      .add_property("raw_pixels_roi",
                     make_getter(&simtbx::nanoBragg::diffBragg::raw_pixels_roi,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::raw_pixels_roi,dcp()),
                    "raw pixels from region of interest only")

      .add_property("oversample_omega",
                     make_getter(&simtbx::nanoBragg::diffBragg::oversample_omega,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::oversample_omega,dcp()),
                    "whether to use an average solid angle correction per pixel, or one at the sub pixel level")

      .add_property("only_save_omega_kahn",
                     make_getter(&simtbx::nanoBragg::diffBragg::only_save_omega_kahn,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::only_save_omega_kahn,dcp()),
                    "ONLY simulate the kahn polarization correction and solid angle components")

      .add_property("max_I_hkl",
                     make_getter(&simtbx::nanoBragg::diffBragg::max_I_hkl,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::max_I_hkl,dcp()),
                    "HKL corresponding to the maximum simulated intensity in the ROI")

      .add_property("Umatrix",
             make_function(&get_U,rbv()),
             make_function(&set_U,dcp()),
             "rotation portion of Amatrix (dxtbx format, C.get_U())")

      .add_property("Bmatrix",
             make_function(&get_B,rbv()),
             make_function(&set_B,dcp()),
             "unit cell portion of Amatrix (dxtbx format; C.get_B())")

      .add_property("Omatrix",
             make_function(&get_O,rbv()),
             make_function(&set_O,dcp()),
             "change of basis Omatrix")

      .add_property("Ncells_abc",
             make_function(&get_Nabc,rbv()),
             make_function(&set_Nabc,dcp()),
             "Get number of mosaic domains along an axis")

      .add_property("update_oversample_during_refinement",
                     make_getter(&simtbx::nanoBragg::diffBragg::update_oversample_during_refinement,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::update_oversample_during_refinement,dcp()),
                     "Allows oversample to change as Ncells abc changes")


    ; // end of diffBragg extention

  } // end of diffBragg_init_module
} // end of namespace
} // end of boost python namespace
} // end of refine
} // end of simtbx

BOOST_PYTHON_MODULE(simtbx_diffBragg_ext)
{
  simtbx::diffBragg::boost_python::diffBragg_init_module();
}
