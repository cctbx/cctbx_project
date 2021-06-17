
#include <scitbx/graphics_utils/colors.h>

#include <boost/python/list.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

namespace scitbx { namespace graphics_utils {

  boost::python::list
  flt_roundoffvec3(
    scitbx::vec3<double> const& vec,
    int const& precision
  )
  {
    // fast version of libtbx.math_utils.roundoff() for a flex vec3 double array. Returns a pythonlist
    boost::python::list retlst;
    retlst.append(flt_roundoff(vec[0], precision));
    retlst.append(flt_roundoff(vec[1], precision));
    retlst.append(flt_roundoff(vec[2], precision));
    return retlst;
  }

namespace {

  void init_module ()
  {
    using namespace boost::python;
    def("flt_roundoff", flt_roundoff, (
      arg("val"),
      arg("precision") = 3
      ));
    def("flt_roundoffvec3", flt_roundoffvec3, (
      arg("val"),
      arg("precision") = 3
      ));
    def("NoNansvec3", NoNansvec3, (
      arg("vecs"),
      arg("defx") = 0.0,
      arg("defy") = 0.0,
      arg("defz") = 0.0
      ));
    def("NoNansHL", NoNansHL, (
      arg("HL"),
      arg("a") = 0.0,
      arg("b") = 0.0,
      arg("c") = 0.0,
      arg("d") = 0.0
      ));
    def("IsNansvec3", IsNansvec3, (
      arg("vecs")));
    def("NoNans", NoNans, (
      arg("arr"),
      arg("def") = 0.0));
    def("IsNans", IsNans, (
      arg("arr")));

    def("make_rainbow_gradient", make_rainbow_gradient, (
      arg("nbins")));
    def("color_rainbow", color_rainbow, (
      arg("selection"),
      arg("color_all")=false));
    def("color_by_property_", color_by_property, (
      arg("properties"),
      arg("selection"),
      arg("color_all") = true,
      arg("gradient_type") = 0,
      arg("min_value") = 0.1));
    def("map_to_rgb_colourmap_", map_to_rgb_colourmap, (
      arg("data_for_colors"),
      arg("colourmap"),
      arg("selection"),
      arg("attenuation"),
      arg("powscale"),
      arg("map_directly"),
      arg("color_all") = true));
    def("color_by_phi_fom_", color_by_phi_fom, (
      arg("phases"),
      arg("foms")));
    def("grayscale_by_property", grayscale_by_property, (
      arg("properties"),
      arg("selection"),
      arg("shade_all")=false,
      arg("invert")=false,
      arg("max_value")=0.95,
      arg("max_value_inverted")=0.1));
    def("scale_selected_colors", scale_selected_colors, (
      arg("input_colors"),
      arg("selection"),
      arg("scale")=0.5));
  }
}
}}

BOOST_PYTHON_MODULE(scitbx_graphics_utils_ext)
{
  scitbx::graphics_utils::init_module();
}
