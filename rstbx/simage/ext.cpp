#include <boost/python.hpp>
#include <rstbx/simage/image_simple.hpp>

namespace rstbx { namespace simage { namespace ext {

  struct image_simple_wrappers
  {
    typedef image_simple wt;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef return_internal_reference<> rir;
      class_<wt>("image_simple", no_init)
        .def(init<bool, bool, bool, bool, bool, bool, bool>((
          arg("apply_proximity_filter")=true,
          arg("apply_detector_clipping")=true,
          arg("apply_proximity_factor")=true,
          arg("store_miller_index_i_seqs")=false,
          arg("store_spots")=false,
          arg("store_signals")=false,
          arg("set_pixels")=false)))
        .def("compute", &wt::compute, (
          arg("unit_cell"),
          arg("miller_indices"),
          arg("spot_intensity_factors"),
          arg("crystal_rotation_matrix"),
          arg("ewald_radius"),
          arg("ewald_proximity"),
          arg("signal_max"),
          arg("detector_distance"),
          arg("detector_size"),
          arg("detector_pixels"),
          arg("point_spread"),
          arg("gaussian_falloff_scale")), rir())
        .add_property("miller_index_i_seqs",
          make_getter(&wt::miller_index_i_seqs, rbv()))
        .add_property("spots", make_getter(&wt::spots, rbv()))
        .add_property("signals", make_getter(&wt::signals, rbv()))
        .def_readonly("pixels", &wt::pixels)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;
    image_simple_wrappers::wrap();
  }

}}} // namespace rstbx::simage::ext

BOOST_PYTHON_MODULE(rstbx_simage_ext)
{
  rstbx::simage::ext::init_module();
}
