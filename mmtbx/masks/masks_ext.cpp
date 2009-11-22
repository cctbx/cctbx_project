#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/def.hpp>
#include <mmtbx/masks/around_atoms.h>
#include "atom_mask.h"
#include "util.h"

namespace mmtbx { namespace masks {
namespace {

  template <typename DataType, typename FloatType>
  struct around_atoms_wrappers
  {
    typedef around_atoms<DataType, FloatType> w_t;

    static
    void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("around_atoms", no_init)
        .def(init<
          cctbx::uctbx::unit_cell const&,
          std::size_t,
          af::shared<scitbx::vec3<double> > const&,
          af::shared<double> const&,
          af::c_grid<3>::index_type const&,
          FloatType const&,
          FloatType const&,
          optional< bool, bool> >((
            arg("unit_cell"),
            arg("space_group_order_z"),
            arg("sites_frac"),
            arg("atom_radii"),
            arg("gridding_n_real"),
            arg("solvent_radius"),
            arg("shrink_truncation_radius"),
            arg("explicit_distance"),
            arg("debug")     )))
        .def_readonly("solvent_radius",
                 &w_t::solvent_radius)
        .def_readonly("n_atom_points",
                 &w_t::n_atom_points)
        .def_readonly("shrink_truncation_radius",
                 &w_t::shrink_truncation_radius)
        .add_property("data", make_getter(&w_t::data, rbv()))
        .def_readonly("contact_surface_fraction",
                 &w_t::contact_surface_fraction)
        .def_readonly("accessible_surface_fraction",
                 &w_t::accessible_surface_fraction)
      ;
    }
  };

  scitbx::af::shared<std::string> generate_groups_p(const std::string &s, int n)
  {
    std::set<std::string> g;
    generate_groups(g, s,n);
    scitbx::af::shared<std::string> result;
    for(std::set<std::string>::const_iterator it=g.begin(); it!=g.end(); ++it)
      result.push_back(*it);
    return result;
  }

  void wrap_atom_mask()
  {
    typedef atom_mask w_t;
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;

    class_<w_t>("atom_mask", no_init)
      .def(init<
          const cctbx::uctbx::unit_cell &,
          const cctbx::sgtbx::space_group &,
          double,
          optional<
            double,
            double,
            double >
          > ((
              arg("unit_cell"),
              arg("group"),
              arg("resolution"),
              arg("grid_step_factor"),
              arg("solvent_radius"),
              arg("shrink_truncation_radius")
              )))
      .def(init<
          const cctbx::uctbx::unit_cell &,
          const cctbx::sgtbx::space_group &,
          const grid_t::index_type &,
          double,
          double
          > ((
              arg("unit_cell"),
              arg("space_group"),
              arg("gridding_n_real"),
              arg("solvent_radius"),
              arg("shrink_truncation_radius")
              )))
      .def("compute", &w_t::compute, (arg("sites_frac"), arg("atom_radii"),
         arg("shells") = shells_array_t() ) )
      .def("structure_factors", &w_t::structure_factors, (arg("indices"),
         arg("layer")=0) )
      .def("grid_size", &w_t::grid_size)
      .def("n_asu_atoms", &w_t::n_asu_atoms)
      .def("n_solvent_layers", &w_t::n_solvent_layers)
      .def("xplor_write_map", &w_t::xplor_write_map, (arg("file_name"),
         arg("layer")=0, arg("invert")=false) )
      .def_readonly("solvent_radius",
               &w_t::solvent_radius)
      .def_readonly("shrink_truncation_radius",
               &w_t::shrink_truncation_radius)
      .def_readonly("contact_surface_fraction",
               &w_t::contact_surface_fraction)
      .def_readonly("accessible_surface_fraction",
               &w_t::accessible_surface_fraction)
      // DO NOTE USE, these are temporarily here
      .def_readonly("debug_mask_asu_time", &w_t::debug_mask_asu_time)
      .def_readonly("debug_atoms_to_asu_time", &w_t::debug_atoms_to_asu_time)
      .def_readonly("debug_accessible_time", &w_t::debug_accessible_time)
      .def_readonly("debug_contact_time", &w_t::debug_contact_time)
      .def_readonly("debug_fft_time", &w_t::debug_fft_time)
      .def_readonly("debug_has_enclosed_box", &w_t::debug_has_enclosed_box)
     ;

    // DO NOT USE
    boost::python::def("generate_groups", generate_groups_p);
  }

  void init_module()
  {
    around_atoms_wrappers<int, double>::wrap();
    wrap_atom_mask();
  }

} // namespace <anonymous>
}} // namespace mmtbx::masks

BOOST_PYTHON_MODULE(mmtbx_masks_ext)
{
  mmtbx::masks::init_module();
}
