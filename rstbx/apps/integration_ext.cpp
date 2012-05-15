#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python.hpp>
#include <rstbx/apps/stills/simple_integration.h>

namespace rstbx { namespace integration { namespace ext {

  struct integration_wrappers
  {

    static void
    wrap()
    {
      using namespace boost::python;

      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;

      class_<simple_integration>(
        "simple_integration", init<>())
         .enable_pickling()
        //Could not figure out how to expose member data from the C++ class into
        // a derived class in Python
        .def("set_pixel_size",&simple_integration::set_pixel_size)
        .def("set_detector_size",&simple_integration::set_detector_size)
        .def("set_frame",&simple_integration::set_frame)
        .def("set_background_factor",&simple_integration::set_background_factor)
        .def("set_nbr_cutoff_sq",&simple_integration::set_nbr_cutoff_sq)
        .def("get_bsmask",&simple_integration::get_bsmask)
        .def("get_ISmask",&simple_integration::get_ISmask)
        .def("positional_correction_mapping",
             &simple_integration::positional_correction_mapping,(
           arg_("predicted"),
           arg_("correction_vectors"),
           arg_("PS_adapt"),
           arg_("IS_adapt"),
           arg_("spots")
            ))
        .def("safe_background",
           (scitbx::af::shared<scitbx::vec2<double> >(simple_integration::*)
           (scitbx::af::shared<scitbx::vec3<double> >,
            annlib_adaptbx::AnnAdaptor const&,
            scitbx::af::shared<int >)
           )
           &simple_integration::safe_background,(
           arg_("predicted"),
           arg_("OS_adapt"),
           arg_("sorted")
            ))
        .def("safe_background",
           (scitbx::af::shared<scitbx::vec2<double> >(simple_integration::*)
           (scitbx::af::shared<scitbx::vec3<double> >,
            annlib_adaptbx::AnnAdaptor const&,
            scitbx::af::shared<int >,
            scitbx::af::shared<int >,
            scitbx::af::shared<int >)
           )
           &simple_integration::safe_background,(
           arg_("predicted"),
           arg_("OS_adapt"),
           arg_("sorted"),
           arg_("tiles"),
           arg_("tile_id")
            ))
        .def("append_ISmask",&simple_integration::append_ISmask)
        .def("integration_proper_fast",
             &simple_integration::integration_proper_fast,(
           arg_("rawdata"),
           arg_("predicted"),
           arg_("hkllist"),
           arg_("detector_xy_draft")
            ))
        .def("get_integrated_data",&simple_integration::get_integrated_data)
        .def("get_integrated_sigma",&simple_integration::get_integrated_sigma)
        .def("get_integrated_miller",&simple_integration::get_integrated_miller)
        .def("get_detector_xy",&simple_integration::get_detector_xy)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;
    integration_wrappers::wrap();
  }

}}} // namespace rstbx::integration::ext

BOOST_PYTHON_MODULE(rstbx_integration_ext)
{
  rstbx::integration::ext::init_module();
}
