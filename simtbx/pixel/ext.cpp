#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <simtbx/pixel/pixel_stats.h>

namespace simtbx { namespace pixel {

  struct pixel_stats_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<simtbx::pixel::pixel_stats>("pixel_stats",init<>() )
        .def(init< >())
        .def("set_whitelist",
             &simtbx::pixel::pixel_stats::set_whitelist,
             (arg_("lunus_filtered_data"), arg_("exp_data"), arg_("sim_mock"))
            )
        .def("set_shoebox_iterator",
             &simtbx::pixel::pixel_stats::set_shoebox_iterator,
             (arg_("shoebox_offset"), arg_("shoebox_size"), arg_("spots_pixels"))
            )
        .def("analyze3",
             &simtbx::pixel::pixel_stats::analyze3,
             (arg_("whitelist_pixels"), arg_("reference_shoebox_sums"),
              arg_("slow_size"), arg_("panel_size"), arg_("keV_per_photon"))
            )
        .def("get_LLG",&simtbx::pixel::pixel_stats::get_LLG)
        .def("get_mnz",&simtbx::pixel::pixel_stats::get_mnz)
        .def("get_sgz",&simtbx::pixel::pixel_stats::get_sgz)
        .def("get_proposal_center_of_mass",&simtbx::pixel::pixel_stats::get_proposal_center_of_mass)
        ;
    }
  };

  } // namespace pixel

  BOOST_PYTHON_MODULE(simtbx_pixel_ext)
  {
    simtbx::pixel::pixel_stats_wrapper::wrap();
  }
} // namespace simtbx
