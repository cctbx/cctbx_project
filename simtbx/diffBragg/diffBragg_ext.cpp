#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python.hpp>
#include <simtbx/diffBragg/diffBragg.h>
//#include <boost/python/tuple.hpp>

using namespace boost::python;
namespace simtbx{
namespace diffBragg{
namespace boost_python { namespace {

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

      .def("get_derivative_pixels", &simtbx::nanoBragg::diffBragg::get_derivative_pixels, "gets the manager raw image")

      .add_property("region_of_interest",
             make_function(&get_roi,rbv()),
             make_function(&set_roi,dcp()),
             "region of interst on detector: fast_min fast_max slow_min slow_max")

      .add_property("raw_pixels_roi",
                     make_getter(&simtbx::nanoBragg::diffBragg::raw_pixels_roi,rbv()),
                     make_setter(&simtbx::nanoBragg::diffBragg::raw_pixels_roi,dcp()),
                     "raw pixels from region of interest only")


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
