#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <iotbx/pdb/hierarchy_v2_atoms.h>

namespace iotbx { namespace pdb { namespace hierarchy_v2 { namespace atoms {

namespace {

  BOOST_PYTHON_FUNCTION_OVERLOADS(reset_tmp_overloads, reset_tmp, 1, 3)

} // namespace <anonymous>

  void
  bpl_wrap()
  {
    using namespace boost::python;
    scitbx::af::boost_python::shared_wrapper<atom>::wrap("af_shared_atom")
      .def("extract_xyz", extract_xyz)
      .def("extract_sigxyz", extract_sigxyz)
      .def("extract_occ", extract_occ)
      .def("extract_sigocc", extract_sigocc)
      .def("extract_b", extract_b)
      .def("extract_sigb", extract_sigb)
      .def("extract_uij", extract_uij)
      .def("extract_siguij", extract_siguij)
      .def("extract_hetero", extract_hetero)
      .def("reset_tmp", reset_tmp, reset_tmp_overloads((
        arg_("self"),
        arg_("first_value")=0,
        arg_("increment")=1)))
      .def("reset_tmp_for_occupancy_groups_simple",
        reset_tmp_for_occupancy_groups_simple)
    ;
  }

}}}} // namespace iotbx::pdb::hierarchy_v2::atoms
