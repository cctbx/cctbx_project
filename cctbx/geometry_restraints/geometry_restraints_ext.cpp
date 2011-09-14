#include <boost/python/module.hpp>

namespace cctbx { namespace geometry_restraints { namespace boost_python {

  void wrap_bond();
  void wrap_bond_similarity();
  void wrap_bond_sorted();
  void wrap_nonbonded();
  void wrap_nonbonded_sorted();
  void wrap_angle();
  void wrap_dihedral();
  void wrap_chirality();
  void wrap_planarity();
  void wrap_motif();

namespace {

  void init_module()
  {
    wrap_bond();
    wrap_bond_similarity();
    wrap_bond_sorted();
    wrap_nonbonded();
    wrap_nonbonded_sorted();
    wrap_angle();
    wrap_dihedral();
    wrap_chirality();
    wrap_planarity();
    wrap_motif();
  }

} // namespace <anonymous>
}}} // namespace cctbx::geometry_restraints::boost_python

BOOST_PYTHON_MODULE(cctbx_geometry_restraints_ext)
{
  cctbx::geometry_restraints::boost_python::init_module();
}
