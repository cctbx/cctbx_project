#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/enum.hpp>

#include <iostream>
#include <rstbx/simulation/sim_xfel1.h>

namespace rstbx { namespace boost_python { namespace {

  void
  init_simulation_module() {
    using namespace boost::python;

    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

    class_<xfel1>("xfel1",init<>())
      .def("set_indices",&xfel1::set_indices)
      .def("set_intensities",&xfel1::set_intensities)
      .def("select_proximal_indices",&xfel1::select_proximal_indices,
        (arg_("half_edge"),
         arg_("detector_distance_m"),
         arg_("pixel_size_m"),
         arg_("orientation"),
         arg_("mosaicity_full_width"),
         arg_("bandpass_full_width"),arg_("wavelength_m"),
         arg_("limiting_resolution_Ang")))
      .def("raw_diffraction",&xfel1::raw_diffraction,
        (arg_("selection"),arg_("pixels"),
         arg_("mosaic_domains"),
         arg_("detector_distance_m"),
         arg_("pixel_size_m"),
         arg_("darwin_factor")
         ))
      .add_property("indices_all",make_getter(&xfel1::indices_all, rbv()))
      .add_property("intensities_all",make_getter(&xfel1::intensities_all, rbv()))
      .add_property("spots", make_getter(&xfel1::spots, rbv()))
      .add_property("signals", make_getter(&xfel1::selection_raw_counts, rbv()))
      .add_property("partialities", make_getter(&xfel1::selection_partiality, rbv()))
    ;
  }

}}} // namespace rstbx::boost_python::<anonymous>


BOOST_PYTHON_MODULE(rstbx_simulation_ext)
{
  rstbx::boost_python::init_simulation_module();
}
