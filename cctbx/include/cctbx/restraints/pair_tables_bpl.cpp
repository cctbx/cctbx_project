#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/stl/map_wrapper.h>
#include <scitbx/stl/vector_wrapper.h>
#include <cctbx/restraints/pair_tables.h>

namespace cctbx { namespace restraints {
namespace {

  struct pair_sym_table_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      scitbx::stl::boost_python::map_wrapper<pair_sym_dict, rir>::wrap(
        "pair_sym_dict");
      scitbx::af::boost_python::shared_wrapper<pair_sym_dict, rir>::wrap(
        "pair_sym_table");
    }
  };

  struct pair_asu_table_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      scitbx::stl::boost_python::map_wrapper<pair_asu_dict, rir>::wrap(
        "pair_asu_dict");
      scitbx::af::boost_python::shared_wrapper<pair_asu_dict, rir>::wrap(
        "pair_asu_table");
    }
  };

  void
  wrap_all()
  {
    pair_sym_table_wrappers::wrap();
    pair_asu_table_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_pair_tables() { wrap_all(); }

}}} // namespace cctbx::restraints::boost_python
