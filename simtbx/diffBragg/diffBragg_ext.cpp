#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python.hpp>
#include <simtbx/diffBragg/diffBragg.h>

using namespace boost::python;
namespace simtbx{
namespace diffBragg{
namespace boost_python { namespace {

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

      .add_property("thetaX",
                 make_getter(&simtbx::nanoBragg::diffBragg::thetaX,rbv()),
                 make_setter(&simtbx::nanoBragg::diffBragg::thetaX,dcp()),
                 "Rotation perturbation about lab X ")

      .add_property("thetaY",
                 make_getter(&simtbx::nanoBragg::diffBragg::thetaY,rbv()),
                 make_setter(&simtbx::nanoBragg::diffBragg::thetaY,dcp()),
                 "Rotation perturbation about lab Y ")

      .add_property("thetaZ",
                 make_getter(&simtbx::nanoBragg::diffBragg::thetaZ,rbv()),
                 make_setter(&simtbx::nanoBragg::diffBragg::thetaZ,dcp()),
                 "Rotation perturbation about lab Z ")


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
