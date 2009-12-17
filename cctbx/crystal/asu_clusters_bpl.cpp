#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_arg.hpp>
#include <cctbx/crystal/asu_clusters.h>
#include <cctbx/crystal/workarounds_bpl.h>

namespace cctbx { namespace crystal {
namespace {

  struct asu_clusters_wrappers
  {
    typedef asu_clusters w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("asu_clusters", no_init)
        .def(init<pair_asu_table<> const&, bool>((
          arg("pair_asu_table"), arg("strictly_in_asu")=true)))
        .def("sort_index_groups_by_size",
          &w_t::sort_index_groups_by_size, return_self<>())
        .def("sort_indices_in_each_group",
          &w_t::sort_indices_in_each_group, return_self<>())
        .add_property("index_groups", make_getter(&w_t::index_groups, rbv()))
      ;
    }
  };

  void
  wrap_all()
  {
    asu_clusters_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_asu_clusters() { wrap_all(); }

}}} // namespace cctbx::crystal::boost_python
