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

  std::string
  combine_rgb_images(
    boost::python::list const& rgb_images)
  {
    namespace bp = boost::python;
    TBXX_ASSERT(bp::len(rgb_images) > 0);
    std::size_t img_size = static_cast<std::size_t>(bp::len(rgb_images[0]));
    std::string result(img_size, '\0');
    std::size_t n_imgs = static_cast<std::size_t>(bp::len(rgb_images));
    boost::scoped_array<const char*> img_ptrs(new const char*[n_imgs]);
    for(std::size_t j=0;j!=n_imgs;j++) {
      TBXX_ASSERT(bp::len(rgb_images[j]) == img_size);
      img_ptrs[j] = bp::extract<const char*>(rgb_images[j])();
    }
    for(std::size_t i=0;i!=img_size;i++) {
      unsigned sum = 0;
      for(std::size_t j=0;j!=n_imgs;j++) {
        sum += static_cast<unsigned>(img_ptrs[j][i]);
      }
      unsigned c = static_cast<unsigned>(
        static_cast<double>(sum) / n_imgs + 0.5);
      if (c > 255) c = 255;
      result[i] = static_cast<char>(c);
    }
    return result;
  }

  void init_module()
  {
    using namespace boost::python;
    image_simple_wrappers::wrap();
    def("combine_rgb_images", combine_rgb_images);
  }

}}} // namespace rstbx::simage::ext

BOOST_PYTHON_MODULE(rstbx_simage_ext)
{
  rstbx::simage::ext::init_module();
}
