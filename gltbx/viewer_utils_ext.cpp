// this module refers throughout to "bonds" and "atoms", but the functions
// actually assume nothing about molecular structures and could be used
// for any similar object.  the one exception is the handling of hydrogens,
// which are often "implicit" in calculations and graphical display,
// whether or not explicit atom objects are present.

#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/array_family/shared.h>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/optional.hpp>

#include <scitbx/vec3.h>

#include <gltbx/error.h>

#include <set>
#include <cmath>

namespace gltbx { namespace viewer_utils {

  namespace af = scitbx::af;

  // Hue, Saturation, Value --> Red, Green, Blue
  scitbx::vec3<double>
  hsv2rgb (double h, double s, double v)
  {
    if (s == 0) {
      return scitbx::vec3<double>(v, v, v);
    }
    h /= 60.0;
    v *= 255.0;
    int i = std::floor(h);
    double f = h - i;
    double p = v * (1 - s);
    double q = v * (1 - (s * f));
    double t = v * (1 - (s * (1 - f)));
    switch (i) {
      case 0  : return scitbx::vec3<double>(v, t, p) / 255.0; break;
      case 1  : return scitbx::vec3<double>(q, v, p) / 255.0; break;
      case 2  : return scitbx::vec3<double>(p, v, t) / 255.0; break;
      case 3  : return scitbx::vec3<double>(p, q, v) / 255.0; break;
      case 4  : return scitbx::vec3<double>(t, p, v) / 255.0; break;
      default : break;
    }
    return scitbx::vec3<double>(v, p, q) / 255.0;
  }

  // this function may be superfluous here, but could be useful elsewhere
  af::shared< scitbx::vec3<double> >
  make_rainbow_gradient (unsigned nbins)
  {
    af::shared< scitbx::vec3<double> > color_gradient(nbins);
    double f_nbins(nbins);
    for (unsigned i = 0; i < nbins; i++) {
      double gradient_ratio = i / f_nbins;
      color_gradient[i] = hsv2rgb(240.0 - (240 * gradient_ratio), 1., 1.);
    }
    return color_gradient;
  }

  // gradient generation is inlined for speed.
  af::shared< scitbx::vec3<double> >
  color_rainbow (
    af::const_ref< bool > const& atoms_visible,
    unsigned visible_atom_count,
    bool color_invisible_atoms=false)
  {
    af::shared< scitbx::vec3<double> > atom_colors(atoms_visible.size());
    unsigned j = 0;
    double f_atom_count = (double) visible_atom_count;
    for (unsigned i_seq = 0; i_seq < atoms_visible.size(); i_seq++) {
      if ((atoms_visible[i_seq]) || (color_invisible_atoms)) {
        double gradient_ratio = j / f_atom_count;
        atom_colors[i_seq] = hsv2rgb(240.0 - (240 * gradient_ratio), 1., 1.);
        j++;
      }
    }
    return atom_colors;
  }

  af::shared< scitbx::vec3<double> >
  scale_selected_colors (
    af::const_ref< scitbx::vec3<double> > const& input_colors,
    af::const_ref< bool > const& selection,
    double scale=0.5)
  {
    GLTBX_ASSERT(input_colors.size() == selection.size());
    GLTBX_ASSERT(scale >= 0);
    af::shared< scitbx::vec3<double> > atom_colors(input_colors.size());
    for (unsigned i_seq = 0; i_seq < input_colors.size(); i_seq++) {
      scitbx::vec3<double> c = input_colors[i_seq];
      if (selection[i_seq]) {
        c[0] *= scale;
        c[1] *= scale;
        c[2] *= scale;
      }
      atom_colors[i_seq] = c;
    }
    return atom_colors;
  }

  af::shared< scitbx::vec3<double> >
  color_by_property (
    af::const_ref< double > const& atom_properties,
    af::const_ref< bool > const& atoms_visible,
    bool color_invisible_atoms=false,
    bool use_rb_color_gradient=false)
  {
    GLTBX_ASSERT(atom_properties.size() > 0);
    af::shared <scitbx::vec3<double> > atom_colors(atom_properties.size());
    double vmax = atom_properties[0];
    double vmin = atom_properties[0];
    for (unsigned i_seq = 0; i_seq < atom_properties.size(); i_seq++) {
      if ((! color_invisible_atoms) && (! atoms_visible[i_seq])) continue;
      if (atom_properties[i_seq] > vmax) vmax = atom_properties[i_seq];
      if (atom_properties[i_seq] < vmin) vmin = atom_properties[i_seq];
    }
    for (unsigned i_seq = 0; i_seq < atom_properties.size(); i_seq++) {
      double gradient_ratio = (atom_properties[i_seq]-vmin) / (vmax-vmin);
      if (use_rb_color_gradient) {
        atom_colors[i_seq] = hsv2rgb(360.0 - (120 * gradient_ratio), 1., 1.);
      } else {
        atom_colors[i_seq] = hsv2rgb(240.0 - (240 * gradient_ratio), 1., 1.);
      }
    }
    return atom_colors;
  }

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

  class atom_visibility
  {
    public :
      af::shared< bool > atoms_visible;
      af::shared< bool > bonds_visible;
      af::shared< bool > points_visible;
      af::shared< bool > selected_bonds_visible;
      af::shared< bool > selected_points_visible;
      unsigned visible_atoms_count;
      unsigned visible_bonds_count;
      unsigned visible_points_count;

      //atom_visibility() {}
      //void get_selection_visibility(
      //  af::const_ref< std::set< unsigned > > const& bonds,
      //  af::const_ref< bool > const& atoms_selected);
      //---end of declarations

      atom_visibility (
        af::const_ref< std::set< unsigned > > const& bonds,
        af::const_ref< bool > atoms_drawable,
        bool flag_show_points)
      {
        typedef std::set< unsigned >::const_iterator it;
        unsigned atom_count = bonds.size();
        visible_atoms_count = 0;
        visible_bonds_count = 0;
        visible_points_count = 0;
        atoms_visible  = af::shared< bool >(atom_count);
        bonds_visible  = af::shared< bool >(atom_count);
        points_visible = af::shared< bool >(atom_count);
        selected_bonds_visible = af::shared< bool >(atom_count); //used later
        selected_points_visible = af::shared< bool >(atom_count); //used later
        for (unsigned i_seq = 0; i_seq < atom_count; i_seq++) {
          if (bonds[i_seq].size() > 0) {
            if (atoms_drawable[i_seq]) {
              atoms_visible[i_seq] = true;
              visible_atoms_count++;
            }
          } else if (flag_show_points == true) {
            atoms_visible[i_seq] = true;
            visible_atoms_count++;
          }
        }
        for (unsigned i_seq = 0; i_seq < atom_count; i_seq++) {
          if (atoms_visible[i_seq]) {
            it j_seqs_end = bonds[i_seq].end();
            for (it j_seq = bonds[i_seq].begin(); j_seq != j_seqs_end; j_seq++){
              if (atoms_visible[*j_seq]) {
                bonds_visible[i_seq] = true;
                visible_bonds_count++;
              }
            }
            if (! bonds_visible[i_seq]) {
              points_visible[i_seq] = true;
              visible_points_count++;
            }
          }
        }
        return;
      }

      void
      get_selection_visibility(
        af::const_ref< std::set< unsigned > > const& bonds,
        af::const_ref< bool > const& atoms_selected)
      {
        typedef std::set< unsigned >::const_iterator it;
        for (unsigned i_seq = 0; i_seq < atoms_selected.size(); i_seq++) {
          if (! atoms_selected[i_seq]) continue;
          if (bonds_visible[i_seq]) {
            it j_seqs_end = bonds[i_seq].end();
            for (it j_seq = bonds[i_seq].begin(); j_seq!=j_seqs_end; j_seq++) {
              if (atoms_selected[*j_seq] && bonds_visible[*j_seq]) {
                selected_bonds_visible[i_seq] = true;
              }
            }
            if (! selected_bonds_visible[i_seq]) {
              selected_points_visible[i_seq] = true;
            }
          } else if (points_visible[i_seq]) {
            selected_points_visible[i_seq] = true;
          }
        }
      }
  };

  scitbx::vec3<double> capped_midpoint (
    scitbx::vec3<double> const& pi,
    scitbx::vec3<double> const& pj,
    double scale_by=1.05) // TODO: scale with view zoom instead
  {
    scitbx::vec3<double> midpoint = (pi + pj) / 2.0;
    scitbx::vec3<double> dxyz = midpoint - pi;
    double length = dxyz.length();
    if (length == 0) {
      return midpoint;
    }
    scitbx::vec3<double> norm_vec = dxyz.normalize();
    scitbx::vec3<double> scaled_vec(norm_vec * length * scale_by);
    return pi + scaled_vec;
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
    typedef std::set< unsigned >::const_iterator it;
    for (unsigned i_seq = 0; i_seq < points.size(); i_seq++) {
      if (! bonds_visible[i_seq]) continue;
      scitbx::vec3<double> const& pi = points[i_seq];
      it j_seqs_end = bonds[i_seq].end();
      for (it j_seq = bonds[i_seq].begin(); j_seq != j_seqs_end; j_seq++) {
        if (! bonds_visible[*j_seq]) continue;
        scitbx::vec3<double> const& pj = points[*j_seq];
        scitbx::vec3<double> const& c = atom_colors[i_seq];
        glColor3f(c[0], c[1], c[2]);
        scitbx::vec3<double> midpoint = capped_midpoint(pi, pj);
        glBegin(GL_LINES);
        glVertex3f(pi[0], pi[1], pi[2]);
        glVertex3f(midpoint[0], midpoint[1], midpoint[2]);
        glEnd();
        //glBegin(GL_POINTS);
        //glVertex3f(midpoint[0], midpoint[1], midpoint[2]);
        //glEnd();
      }
    }
    handle_error();
  }

  double line_given_points_distance_sq (
    scitbx::vec3<double> point,
    scitbx::vec3<double> reference_point,
    scitbx::vec3<double> delta,
    double delta_norm_sq)
  {
    if (delta_norm_sq == 0) {
      return (point - reference_point).length_sq();
    }
    return delta.cross(point - reference_point).length_sq() / delta_norm_sq;
  }

  boost::optional<unsigned>
  closest_visible_point (
    af::const_ref< scitbx::vec3<double> > const& points,
    af::const_ref< bool > const& atoms_visible,
    scitbx::vec3<double> const& point0,
    scitbx::vec3<double> const& point1,
    double min_dist_sq = 1)
  {
    GLTBX_ASSERT(atoms_visible.size() == points.size());
    scitbx::vec3<double> delta = point1 - point0;
    double delta_norm_sq = delta.length_sq();
    bool found_neighbor_point = false;
    unsigned closest_i_seq = 0;
    for (unsigned i_seq = 0; i_seq < points.size(); i_seq++) {
      if (! atoms_visible[i_seq]) continue;
      double dist_sq = line_given_points_distance_sq(points[i_seq],
        point0, delta, delta_norm_sq);
      if (min_dist_sq > dist_sq) {
        min_dist_sq = dist_sq;
        closest_i_seq = i_seq;
        found_neighbor_point = true;
      }
    }
    if (found_neighbor_point) {
      return boost::optional<unsigned>(closest_i_seq);
    }
    return boost::optional<unsigned>();
  }

  void
  init_module()
  {
    using namespace boost::python;
    def("make_rainbow_gradient", make_rainbow_gradient, (
      arg("nbins")));
    def("color_rainbow", color_rainbow, (
      arg("atoms_visible"),
      arg("visible_atom_count"),
      arg("color_invisible_atoms")=false));
    def("color_by_property", color_by_property, (
      arg("atom_properties"),
      arg("atoms_visible"),
      arg("color_invisible_atoms")=true,
      arg("use_rb_color_gradient")=false));
    def("scale_selected_colors", scale_selected_colors, (
      arg("input_colors"),
      arg("selection"),
      arg("scale")=0.5));
    def("draw_points", draw_points, (
      arg("points"),
      arg("atom_colors"),
      arg("points_visible"),
      arg("cross_radius")=0.25));
    def("draw_bonds", draw_bonds, (
      arg("points"),
      arg("bonds"),
      arg("atom_colors"),
      arg("bonds_visible")));
    def("closest_visible_point", closest_visible_point, (
      arg("points"),
      arg("atoms_visible"),
      arg("point0"),
      arg("point1"),
      arg("min_dist_sq")=1.0));
    typedef atom_visibility a_v;
    typedef return_value_policy<return_by_value> rbv;
    class_<a_v>("atom_visibility", no_init)
      .def(init<af::const_ref< std::set<unsigned> > const&,
                af::const_ref<bool> const&,
                bool>((
        arg("bonds"),
        arg("atoms_drawable"),
        arg("flag_show_points"))))
      .def("get_selection_visibility", &a_v::get_selection_visibility, (
        arg("bonds"),
        arg("atoms_selected")))
      .def_readonly("visible_atoms_count", &a_v::visible_atoms_count)
      .def_readonly("visible_bonds_count", &a_v::visible_bonds_count)
      .def_readonly("visible_points_count", &a_v::visible_points_count)
      .add_property("atoms_visible", make_getter(&a_v::atoms_visible, rbv()))
      .add_property("bonds_visible", make_getter(&a_v::bonds_visible, rbv()))
      .add_property("points_visible", make_getter(&a_v::points_visible, rbv()))
      .add_property("selected_bonds_visible",
        make_getter(&a_v::selected_bonds_visible, rbv()))
      .add_property("selected_points_visible",
        make_getter(&a_v::selected_points_visible, rbv()))
    ;
  }

}} // namespace gltbx::viewer_utils

BOOST_PYTHON_MODULE(gltbx_viewer_utils_ext)
{
  gltbx::viewer_utils::init_module();
}
