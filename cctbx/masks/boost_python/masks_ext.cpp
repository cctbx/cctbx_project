#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/def.hpp>
#include <cctbx/masks/around_atoms.h>
#include <cctbx/masks/flood_fill.h>

namespace cctbx { namespace masks {
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

  template <typename DataType, typename FloatType>
  struct flood_fill_wrappers
  {
    typedef flood_fill<DataType, FloatType> w_t;

    static
    void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("flood_fill", no_init)
        .def(init<
          af::ref<DataType, af::c_grid_periodic<3> > const & >((arg("data"))))
        .def("n_voids", &w_t::n_voids)
        .def("averaged_indices", &w_t::averaged_indices)
        .def("averaged_frac_coords", &w_t::averaged_frac_coords)
        .def("grid_points_per_void", make_getter(&w_t::grid_points_per_void, rbv()))
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    around_atoms_wrappers<int, double>::wrap();
    flood_fill_wrappers<int, double>::wrap();
  }

} // namespace <anonymous>
}} // namespace cctbx::masks

BOOST_PYTHON_MODULE(cctbx_masks_ext)
{
  cctbx::masks::init_module();
}
