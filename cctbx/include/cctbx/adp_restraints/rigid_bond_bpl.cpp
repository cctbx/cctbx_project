#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/adp_restraints/rigid_bond.h>

namespace cctbx { namespace adp_restraints {
namespace {

  struct rigid_bond_wrappers
  {

    static void
    wrap()
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
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    rigid_bond_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_rigid_bond() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
