#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <cctbx/adp_restraints/rigid_bond.h>

namespace cctbx { namespace adp_restraints {
namespace {

  void init_module()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_;
    class_<rigid_bond_pair>("rigid_bond_pair",
                             init<vec3<double> const&,
                                  vec3<double> const&,
                                  sym_mat3<double> const&,
                                  sym_mat3<double> const&,
                                  cctbx::uctbx::unit_cell const&>())
      .def("z_12", &rigid_bond_pair::z_12)
      .def("z_21", &rigid_bond_pair::z_21)
      .def("delta_z", &rigid_bond_pair::delta_z)
    ;
  }

}
}} // namespace cctbx::adp_restraints

BOOST_PYTHON_MODULE(cctbx_rigid_bond_ext)
{
  cctbx::adp_restraints::init_module();
}
