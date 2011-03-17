
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

  void wrap_simple_restraints()
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

  void wrap_lennard_jones() {
    using namespace boost::python;
    typedef h_bond_lj_proxy w_t;
    class_<w_t>("h_bond_lj_proxy", no_init)
      .def(init<
        af::tiny<unsigned, 2> const&, double, double >((
          arg("i_seqs"),
          arg("distance_ideal"),
          arg("distance_cut"))))
      //.def_readonly("i_seqs", &w_t::i_seqs)
      .def_readonly("distance_ideal", &w_t::distance_ideal)
      .def_readonly("distance_cut", &w_t::distance_cut);
    {
      typedef return_internal_reference<> rir;
      scitbx::af::boost_python::shared_wrapper<h_bond_lj_proxy, rir>::wrap(
        "shared_h_bond_lennard_jones_proxy");
    }

    def("h_bond_lennard_jones_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<h_bond_lj_proxy> const&,
        af::ref<scitbx::vec3<double> > const&,
        double,
        double,
        int,
        int,
        double,
        bool))
      h_bond_lennard_jones_residual_sum, (
      arg("sites_cart"),
      arg("proxies"),
      arg("gradient_array"),
      arg("falloff_distance")=0.05,
      arg("sigma_base")=0.81649658092772603,
      arg("a")=6,
      arg("b")=4,
      arg("scale")=1.0,
      arg("use_finite_differences")=true));
  }

  void wrap_implicit_restraints ()
  {
    using namespace boost::python;
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

  void wrap_misc () {
    using namespace boost::python;
    def("simple_hbonds_as_simple_bonds", simple_hbonds_as_simple_bonds, (
      arg("proxies")));
    def("lj_hbonds_as_simple_bonds", lj_hbonds_as_simple_bonds, (
      arg("proxies")));
    def("implicit_hbonds_as_simple_bonds", implicit_hbonds_as_simple_bonds, (
      arg("proxies")));
  }

} // namespace anonymous

namespace boost_python {

  void wrap_hbond() {
    wrap_simple_restraints();
    wrap_lennard_jones();
    wrap_implicit_restraints();
    wrap_misc();
  }

} // namespace boost_python

}} //namespace mmtbx::geometry_restraints

BOOST_PYTHON_MODULE(mmtbx_hbond_restraints_ext)
{
  mmtbx::geometry_restraints::boost_python::wrap_hbond();
}
