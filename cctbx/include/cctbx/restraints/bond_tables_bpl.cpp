#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/restraints/bond_tables.h>

namespace cctbx { namespace restraints {
namespace {

  struct bond_sym_dict_wrappers
  {
    typedef bond_sym_dict w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t>("bond_sym_dict")
        .def(map_indexing_suite<w_t, true>())
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<bond_sym_dict, rir>::wrap(
          "bond_sym_table");
      }
    }
  };

  void
  wrap_all()
  {
    bond_sym_dict_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_bond_tables() { wrap_all(); }

}}} // namespace cctbx::restraints::boost_python
