
#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/optional.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

#include <mmtbx/geometry_restraints/hbond.h>

namespace mmtbx { namespace geometry_restraints {
namespace {

  void wrap_hbond_restraints()
  {
    using namespace boost::python;
    typedef h_bond_simple_proxy w_t;
    class_<w_t>("h_bond_simple_proxy", no_init)
      .def(init<
        af::tiny<unsigned, 2> const&, double, double, double, double >((
          arg("i_seqs"),
          arg("distance_ideal"),
          arg("distance_cut"),
          arg("weight"),
          arg("slack")=0)))
      //.def_readonly("i_seqs", &w_t::i_seqs)
      .def_readonly("distance_ideal", &w_t::distance_ideal)
      .def_readonly("distance_cut", &w_t::distance_cut)
      .def_readonly("weight", &w_t::weight)
      .def_readonly("slack", &w_t::slack)
    ;
    {
      typedef return_internal_reference<> rir;
      scitbx::af::boost_python::shared_wrapper<h_bond_simple_proxy, rir>::wrap(
        "shared_h_bond_simple_proxy");
    }

    def("h_bond_simple_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<h_bond_simple_proxy> const&,
        af::ref<scitbx::vec3<double> > const&,
        double,
        double))
      h_bond_simple_residual_sum, (
      arg("sites_cart"),
      arg("proxies"),
      arg("gradient_array"),
      arg("hbond_weight")=1.0,
      arg("falloff_distance")=0.05));

  }
} // namespace anonymous

namespace boost_python {

  void wrap_hbond() { wrap_hbond_restraints(); }

} // namespace boost_python

}} //namespace mmtbx::geometry_restraints

BOOST_PYTHON_MODULE(mmtbx_hbond_restraints_ext)
{
  mmtbx::geometry_restraints::boost_python::wrap_hbond();
}
