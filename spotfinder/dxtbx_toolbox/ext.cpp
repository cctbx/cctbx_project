#include <boost/python.hpp>
#include <boost/python/args.hpp>
#include <spotfinder/dxtbx_toolbox/distl.h>

using namespace boost::python;
using namespace spotfinder::dxtbx;
namespace af = scitbx::af;

BOOST_PYTHON_MODULE(spotfinder_dxtbx_ext)
{
   typedef return_value_policy<return_by_value> rbv;
   typedef default_call_policies dcp;


   class_<w_Distl>("w_Distl", no_init)
     .def(init<std::string, bool >((arg_("optionstring"),arg_("report_overloads"))))
     .def("set_resolution_outer",&w_Distl::set_resolution_outer)
     .def("setspotimg",&w_Distl::setspotimg, (
          arg_("panel"),arg_("beam"),
          arg_("rawdata"),
          arg_("peripheral_margin"),arg_("saturation")
      ))
     .def("set_tiling",(void(w_Distl::*)(const string&))&w_Distl::set_tiling)
     .def("set_tiling",(void(w_Distl::*)(af::flex_int const&,int const&))&w_Distl::set_tiling,
         (arg_("detector_tiling"),arg_("peripheral_margin"))
      )
     .def("Z_data",&w_Distl::Z_data)
     .def("mod_data",&w_Distl::mod_data)
     .def("get_underload",&w_Distl::get_underload)
     .def("set_minimum_spot_area",&w_Distl::set_minimum_spot_area)
     .def("get_minimum_spot_area",&w_Distl::get_minimum_spot_area)
     .def("set_minimum_signal_height",&w_Distl::set_minimum_signal_height)
     .def("set_minimum_spot_height",&w_Distl::set_minimum_spot_height)
     .def("set_spot_area_maximum_factor",&w_Distl::set_spot_area_maximum_factor)
     .def("set_peak_intensity_maximum_factor",&w_Distl::set_peak_intensity_maximum_factor)
     .def("set_scanbox_windows",&w_Distl::set_scanbox_windows)
     .def("parameter_guarantees",&w_Distl::parameter_guarantees)
     .def("pxlclassify",&w_Distl::pxlclassify)
     .def("search_icerings",&w_Distl::search_icerings)
     .def("search_maximas",&w_Distl::search_maximas)
     .def("search_spots",&w_Distl::search_spots)
     .def("search_overloadpatches",&w_Distl::search_overloadpatches)
     .def("finish_analysis",&w_Distl::finish_analysis)
     .add_property("spots",make_getter(&w_Distl::spots,rbv()),
                           make_setter(&w_Distl::spots,dcp()))
     .def("nicerings",&w_Distl::nicerings)
     .def("isIsolated",&w_Distl::isIsolated)
     .add_property("icerings",make_getter(&w_Distl::icerings,rbv()))
     .def("imgresol",&w_Distl::imgresol)
     .def("background_resolutions",&w_Distl::background_resolutions)
     .def("background_means",&w_Distl::background_means)
     .def("background_wndw_sz",&w_Distl::background_wndw_sz)
     .def("spotbasesize",&w_Distl::get_minimum_spot_area)
   ;

}
