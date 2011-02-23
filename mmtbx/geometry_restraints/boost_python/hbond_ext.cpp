
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
  using namespace boost::python;

  void wrap_simple_restraints()
  {
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

  void wrap_implicit_restraints ()
  {
    typedef h_bond_implicit_proxy w_t;
    class_<w_t>("h_bond_implicit_proxy", no_init)
      .def(init<
        af::tiny<unsigned,3> const&, double, double, double, double, double >((
          arg("i_seqs"),
          arg("distance_ideal"),
          arg("distance_cut"),
          arg("theta_high"),
          arg("theta_low"),
          arg("weight"))))
      .def_readonly("distance_ideal", &w_t::distance_ideal)
      .def_readonly("distance_cut", &w_t::distance_cut)
      .def_readonly("theta_high", &w_t::theta_high)
      .def_readonly("theta_low", &w_t::theta_low)
      .def_readonly("weight", &w_t::weight)
    ;
    {
      typedef return_internal_reference<> rir;
      scitbx::af::boost_python::shared_wrapper<h_bond_implicit_proxy,rir>::wrap(
        "shared_h_bond_implicit_proxy");
    }

    def("h_bond_implicit_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<h_bond_implicit_proxy> const&,
        af::ref<scitbx::vec3<double> > const&,
        double,
        double,
        bool))
      h_bond_implicit_residual_sum, (
      arg("sites_cart"),
      arg("proxies"),
      arg("gradient_array"),
      arg("falloff_distance")=0.05,
      arg("epsilon")=0.0001,
      arg("use_finite_differences")=true));
  }

} // namespace anonymous

namespace boost_python {

  void wrap_hbond() {
    wrap_simple_restraints();
    wrap_implicit_restraints();
  }

} // namespace boost_python

}} //namespace mmtbx::geometry_restraints

BOOST_PYTHON_MODULE(mmtbx_hbond_restraints_ext)
{
  mmtbx::geometry_restraints::boost_python::wrap_hbond();
}
