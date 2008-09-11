#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/optional.hpp>

#include <scitbx/vec3.h>

#include <gltbx/error.h>

#include <set>

namespace gltbx { namespace viewer_utils {

  namespace af = scitbx::af;


  void
  draw_points (
    af::const_ref< scitbx::vec3<double> > const& points,
    af::const_ref< scitbx::vec3<double> > const& atom_colors,
    af::const_ref< bool > const& points_visible,
    double cross_radius=0.25)
  {
    GLTBX_ASSERT(atom_colors.size() == points.size());
    GLTBX_ASSERT(points_visible.size() == points.size());
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
    handle_error();
  }

  void
  draw_bonds (
    af::const_ref< scitbx::vec3<double> > const& points,
    af::const_ref< std::set< unsigned > > const& bonds,
    af::const_ref< scitbx::vec3<double> > const& atom_colors,
    af::const_ref< bool > const& bonds_visible)
  {
    GLTBX_ASSERT(bonds.size() == points.size());
    GLTBX_ASSERT(atom_colors.size() == points.size());
    GLTBX_ASSERT(bonds_visible.size() == points.size());
    for (unsigned i_seq = 0; i_seq < points.size(); i_seq++) {
      if (! bonds_visible[i_seq]) continue;
      scitbx::vec3<double> const& pi = points[i_seq];
      typedef std::set< unsigned >::const_iterator it;
      it j_seqs_end = bonds[i_seq].end();
      for (it j_seq = bonds[i_seq].begin(); j_seq != j_seqs_end; j_seq++) {
        if (! bonds_visible[*j_seq]) continue;
        scitbx::vec3<double> const& pj = points[*j_seq];
        scitbx::vec3<double> const& c = atom_colors[i_seq];
        glColor3f(c[0], c[1], c[2]);
        scitbx::vec3<double> midpoint = (pi + pj) / 2.0;
        glBegin(GL_LINES);
        glVertex3f(pi[0], pi[1], pi[2]);
        glVertex3f(midpoint[0], midpoint[1], midpoint[2]);
        glEnd();
      }
    }
    handle_error();
  }

  double norm_sq (
    scitbx::vec3<double> point)
  {
    return (point[0]*point[0])+(point[1]*point[1])+(point[2]*point[2]);
  }

  scitbx::vec3<double> cross (
    scitbx::vec3<double> a,
    scitbx::vec3<double> b)
  {
    scitbx::vec3<double> product;
    product[0] = a[1] * b[1] - b[1] * a[2];
    product[1] = a[2] * b[0] - b[2] * a[0];
    product[2] = a[0] * b[1] - b[0] * a[1];
    return product;
  }

  double distance_sq (
    scitbx::vec3<double> point,
    scitbx::vec3<double> reference_point,
    scitbx::vec3<double> delta,
    double delta_norm_sq)
  {
    if (delta_norm_sq == 0) {
      return norm_sq(point - reference_point);
    }
    scitbx::vec3<double> cross_product = cross(delta, point - reference_point);
    return norm_sq(cross_product) / delta_norm_sq;
  }

  boost::optional<unsigned> closest_visible_point (
    af::const_ref< scitbx::vec3<double> > const& points,
    af::const_ref< bool > const& atoms_visible,
    af::const_ref< scitbx::vec3<double> > const& pick_points,
    double min_dist_sq = 1)
  {
    GLTBX_ASSERT(atoms_visible.size() == points.size());
    GLTBX_ASSERT(pick_points.size() == 2);
    scitbx::vec3<double> reference_point = pick_points[0];
    scitbx::vec3<double> delta = pick_points[1] - reference_point;
    double delta_norm_sq = norm_sq(delta);
    //double min_dist_sq = 1;
    bool found_neighbor_point = false;
    unsigned closest_i_seq = 0;
    for (unsigned i_seq = 0; i_seq < points.size(); i_seq++) {
      if (! atoms_visible[i_seq]) continue;
      double dist_sq = distance_sq(points[i_seq], reference_point, delta,
        delta_norm_sq);
      if (min_dist_sq > dist_sq) {
        min_dist_sq = dist_sq;
        closest_i_seq = i_seq;
        found_neighbor_point = true;
      }
    }
    if (found_neighbor_point) {
      return boost::optional<unsigned>(closest_i_seq);
    } else {
      return boost::optional<unsigned>();
    }
  }

  BOOST_PYTHON_FUNCTION_OVERLOADS(draw_points_overloads, draw_points, 3, 4)
  BOOST_PYTHON_FUNCTION_OVERLOADS(closest_visible_point_overloads,
    closest_visible_point, 3, 4)

  void
  init_module()
  {
    using namespace boost::python;
    def("draw_points", draw_points, draw_points_overloads((
      arg_("points"),
      arg_("atom_colors"),
      arg_("points_visible"),
      arg_("cross_radius")=0.25)));
    def("draw_bonds", draw_bonds, (
      arg_("points"),
      arg_("bonds"),
      arg_("atom_colors"),
      arg_("bonds_visible")));
    def("closest_visible_point", closest_visible_point,
      closest_visible_point_overloads((
        arg_("points"),
        arg_("atoms_visible"),
        arg_("pick_points"),
        arg_("min_dist_sq")=1.0)));
  }

}} // namespace gltbx::viewer_utils

BOOST_PYTHON_MODULE(gltbx_viewer_utils_ext)
{
  gltbx::viewer_utils::init_module();
}
