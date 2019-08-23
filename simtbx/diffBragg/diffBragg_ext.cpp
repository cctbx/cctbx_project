#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python.hpp>
#include <simtbx/diffBragg/diffBragg.h>

using namespace boost::python;
namespace simtbx{
namespace diffBragg{
namespace boost_python { namespace {

  /* spindle starting angle phi, in deg */
  static double get_thetaZ(simtbx::nanoBragg::diffBragg const& diffBragg) {
      return diffBragg.thetaZ*simtbx::nanoBragg::RTD;
  }
  static void   set_thetaZ(simtbx::nanoBragg::diffBragg& diffBragg, double const& value) {
      diffBragg.thetaZ = value/simtbx::nanoBragg::RTD;
  }

  static double get_thetaY(simtbx::nanoBragg::diffBragg const& diffBragg) {
      return diffBragg.thetaY*simtbx::nanoBragg::RTD;
  }
  static void   set_thetaY(simtbx::nanoBragg::diffBragg& diffBragg, double const& value) {
      diffBragg.thetaY = value/simtbx::nanoBragg::RTD;
  }

  static double get_thetaX(simtbx::nanoBragg::diffBragg const& diffBragg) {
      return diffBragg.thetaX*simtbx::nanoBragg::RTD;
  }
  static void   set_thetaX(simtbx::nanoBragg::diffBragg& diffBragg, double const& value) {
      diffBragg.thetaX = value/simtbx::nanoBragg::RTD;
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

      .add_property("thetaZ",
                 make_function(&get_thetaZ,rbv()),
                 make_function(&set_thetaZ,dcp()),
                 "Rotation perturbation about lab Z ")

      .add_property("thetaY",
                 make_function(&get_thetaY,rbv()),
                 make_function(&set_thetaY,dcp()),
                 "Rotation perturbation about lab Y ")

      .add_property("thetaX",
                 make_function(&get_thetaX,rbv()),
                 make_function(&set_thetaX,dcp()),
                 "Rotation perturbation about lab X ")


    ; // end of class def
  } // end of diffBragg_init_module
} // end of namespace
} // end of boost python namespace
} // end of refine
} // end of simtbx

BOOST_PYTHON_MODULE(simtbx_diffBragg_ext)
{
  simtbx::diffBragg::boost_python::diffBragg_init_module();
}
