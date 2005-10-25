#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <mmtbx/masks/around_atoms.h>

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
          FloatType const&>((
            arg_("unit_cell"),
            arg_("space_group_order_z"),
            arg_("sites_frac"),
            arg_("atom_radii"),
            arg_("gridding_n_real"),
            arg_("solvent_radius"),
            arg_("shrink_truncation_radius"))))
        .def_readonly("solvent_radius",
                 &w_t::solvent_radius)
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

  void init_module()
  {
    around_atoms_wrappers<int, float>::wrap();
  }

} // namespace <anonymous>
}} // namespace mmtbx::masks

BOOST_PYTHON_MODULE(mmtbx_masks_ext)
{
  mmtbx::masks::init_module();
}
