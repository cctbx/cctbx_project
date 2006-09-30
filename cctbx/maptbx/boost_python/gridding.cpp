#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/gridding.h>
#include <boost/python/def.hpp>

namespace cctbx { namespace maptbx { namespace boost_python {

  void wrap_gridding()
  {
    using namespace boost::python;

    def("determine_gridding",
      (af::tiny<int, 3>(*)
        (uctbx::unit_cell const&,
         double,
         double,
         af::tiny<int, 3> const&,
         int,
         bool)) determine_gridding);

    def("determine_gridding",
      (af::tiny<int, 3>(*)
        (uctbx::unit_cell const&,
         double,
         double,
         sgtbx::search_symmetry_flags const&,
         sgtbx::space_group_type const&,
         af::tiny<int, 3> const&,
         int,
         bool)) determine_gridding);
  }

}}} // namespace cctbx::maptbx::boost_python
