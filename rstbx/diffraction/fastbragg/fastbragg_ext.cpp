#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/enum.hpp>

#include <rstbx/diffraction/fastbragg/fastbragg.h>

using namespace boost::python;

namespace rstbx {
namespace diffraction {
namespace fastbragg {

namespace boost_python { namespace {

  boost::python::tuple
  foo()
  {
    return boost::python::make_tuple(1,2,3,4);
  }

  void
  init_module() {
    using namespace boost::python;

    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

    def("foo", &foo);

    class_<detector>("detector",no_init)
      .enable_pickling()
      .def(init<int const&, int const&, double const&>(
        (arg_("slowpixels"), arg_("fastpixels"), arg_("pixel_size"))))

      .def("set_region_of_interest",&detector::set_region_of_interest,
        (arg_("slow_min"), arg_("slow_max"),
         arg_("fast_min"), arg_("fast_max")))

      .def("set_oversampling",&detector::set_oversampling,
        (arg_("oversampling")))

      .def("intdata",&detector::intdata)

      .add_property("raw",
                     make_getter(&detector::raw,rbv()))
      .add_property("pixel_sz",
                     make_getter(&detector::pixel_sz,rbv()))
    ;

    class_<camera>("camera",init<>())
      .add_property("distance",
                     make_getter(&camera::distance,rbv()),
                     make_setter(&camera::distance,dcp()))
      .add_property("Ybeam",
                     make_getter(&camera::Ybeam,rbv()),
                     make_setter(&camera::Ybeam,dcp()))
      .add_property("Zbeam",
                     make_getter(&camera::Zbeam,rbv()),
                     make_setter(&camera::Zbeam,dcp()))
      .add_property("lambda0",
                     make_getter(&camera::lambda0,rbv()),
                     make_setter(&camera::lambda0,dcp()))
      .add_property("dispersion",
                     make_getter(&camera::dispersion,rbv()),
                     make_setter(&camera::dispersion,dcp()))
      .add_property("dispsteps",
                     make_getter(&camera::dispsteps,rbv()),
                     make_setter(&camera::dispsteps,dcp()))
      .add_property("hdivrange",
                     make_getter(&camera::hdivrange,rbv()),
                     make_setter(&camera::hdivrange,dcp()))
      .add_property("vdivrange",
                     make_getter(&camera::vdivrange,rbv()),
                     make_setter(&camera::vdivrange,dcp()))
      .add_property("hdivstep",
                     make_getter(&camera::hdivstep,rbv()),
                     make_setter(&camera::hdivstep,dcp()))
      .add_property("vdivstep",
                     make_getter(&camera::vdivstep,rbv()),
                     make_setter(&camera::vdivstep,dcp()))
      .add_property("source_distance",
                     make_getter(&camera::source_distance,rbv()),
                     make_setter(&camera::source_distance,dcp()))
      .add_property("fluence",
                     make_getter(&camera::fluence,rbv()),
                     make_setter(&camera::fluence,dcp()))
      .def("get_wavelengths",&camera::get_wavelengths)
    ;

    class_<crystal>("crystal",init<>())
      .add_property("orientation",
                     make_getter(&crystal::orientation,rbv()),
                     make_setter(&crystal::orientation,dcp()))
      .add_property("miller",
                     make_getter(&crystal::miller,rbv()),
                     make_setter(&crystal::miller,dcp()))
      .add_property("amplitudes",
                     make_getter(&crystal::amplitudes,rbv()),
                     make_setter(&crystal::amplitudes,dcp()))
      .add_property("Na",
                     make_getter(&crystal::Na,rbv()),
                     make_setter(&crystal::Na,dcp()))
      .add_property("Nb",
                     make_getter(&crystal::Nb,rbv()),
                     make_setter(&crystal::Nb,dcp()))
      .add_property("Nc",
                     make_getter(&crystal::Nc,rbv()),
                     make_setter(&crystal::Nc,dcp()))
    ;

    class_<fast_bragg_simulation>("fast_bragg_simulation",init<>())

      .def("set_detector",&fast_bragg_simulation::set_detector)
      .def("set_camera",&fast_bragg_simulation::set_camera)
      .def("set_crystal",&fast_bragg_simulation::set_crystal)
      .def("sweep_over_detector",&fast_bragg_simulation::sweep_over_detector,
        (arg_("verbose")=false))
      .def("to_smv_format",&fast_bragg_simulation::to_smv_format,
        (arg_("fileout"),arg_("intfile_scale")=0,
         arg_("saturation")=65535,arg_("verbose")=false))
    ;
  }

}}}}} // namespace rstbx::diffraction::fastbragg::boost_python::<anonymous>

BOOST_PYTHON_MODULE(rstbx_diffraction_fastbragg_ext)
{
  rstbx::diffraction::fastbragg::boost_python::init_module();
}
