#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/stl/vector_wrapper.h>
#include <cctbx/restraints/bond_tables.h>

namespace cctbx { namespace restraints {
namespace {

  struct bond_params_table_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<bond_params_dict>("bond_params_dict")
        .def(map_indexing_suite<bond_params_dict>())
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<bond_params_dict, rir>::wrap(
          "bond_params_table");
      }
    }
  };

  struct bond_sym_table_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<bond_sym_dict>("bond_sym_dict")
        .def(map_indexing_suite<bond_sym_dict>())
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<bond_sym_dict, rir>::wrap(
          "bond_sym_table");
      }
    }
  };

  struct bond_asu_table_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<bond_asu_dict>("bond_asu_dict")
        .def(map_indexing_suite<bond_asu_dict>())
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<bond_asu_dict, rir>::wrap(
          "bond_asu_table");
      }
    }
  };

  void
  wrap_all()
  {
    bond_params_table_wrappers::wrap();
    bond_sym_table_wrappers::wrap();
    bond_asu_table_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_bond_tables() { wrap_all(); }

}}} // namespace cctbx::restraints::boost_python
