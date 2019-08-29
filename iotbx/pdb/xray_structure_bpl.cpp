#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <iotbx/pdb/xray_structure.h>

namespace iotbx { namespace pdb {
namespace {

  struct xray_structures_simple_extension_wrappers
  {
    typedef xray_structures_simple_extension<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("xray_structures_simple_extension", no_init)
        .def(init<
               bool,
               bool,
               bool,
               bool,
               bool,
               af::shared<hierarchy::atom_with_labels> const&,
               af::shared<std::size_t> const&,
               std::set<std::string> const&,
               cctbx::uctbx::unit_cell const&,
               scitbx::mat3<double> const&,
               scitbx::vec3<double> const&>())
        .add_property("scatterers", make_getter(&w_t::scatterers, rbv()))
        .def("next", &w_t::next)
        .def("__next__", &w_t::next)
      ;
    }
  };

  void
  wrap_xray_structure_impl()
  {
    xray_structures_simple_extension_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_xray_structure() { wrap_xray_structure_impl(); }

}}} // namespace iotbx::pdb::boost_python
