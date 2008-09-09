#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

#include <scitbx/vec3.h>
#include <gltbx/include_opengl.h>

namespace gltbx { namespace viewer {

  namespace af = scitbx::af;


  void
  draw_points (
    af::const_ref< scitbx::vec3<double> > const& points,
    af::const_ref< scitbx::vec3<double> > const& atom_colors,
    af::const_ref< bool > const& points_visible,
    double cross_radius=0.25)
  {
    double f = cross_radius;
    for (unsigned i_seq = 0; i_seq < points.size(); i_seq++) {
      if (! points_visible[i_seq]) continue;
      double x = points[i_seq][0];
      double y = points[i_seq][1];
      double z = points[i_seq][2];
      glBegin(GL_LINES);
      scitbx::vec3<double> const& c = atom_colors[i_seq];
      glColor3f(c[0], c[1], c[2]);
      glVertex3f(x-f,y,z);
      glVertex3f(x+f,y,z);
      glVertex3f(x,y-f,z);
      glVertex3f(x,y+f,z);
      glVertex3f(x,y,z-f);
      glVertex3f(x,y,z+f);
      glEnd();
    }
  }

  BOOST_PYTHON_FUNCTION_OVERLOADS(draw_points_overloads, draw_points, 3, 4)

  void
  init_module()
  {
    using namespace boost::python;
    def("draw_points", draw_points, draw_points_overloads((
      arg_("points"),
      arg_("atom_colors"),
      arg_("points_visible"),
      arg_("cross_radius")=0.25)));
  }

}} // namespace gltbx::viewer

BOOST_PYTHON_MODULE(gltbx_viewer_ext)
{
  gltbx::viewer::init_module();
}
