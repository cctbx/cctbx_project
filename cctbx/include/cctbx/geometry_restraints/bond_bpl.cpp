#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/stl/map_wrapper.h>
#include <cctbx/geometry_restraints/bond.h>

namespace cctbx { namespace geometry_restraints {
namespace {

  struct bond_params_wrappers
  {
    typedef bond_params w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("bond_params", no_init)
        .def(init<double, double>((arg_("distance_ideal"), arg_("weight"))))
        .def_readwrite("distance_ideal", &w_t::distance_ideal)
        .def_readwrite("weight", &w_t::weight)
      ;
    }
  };

  struct bond_params_table_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      scitbx::stl::boost_python::map_wrapper<bond_params_dict, rir>::wrap(
        "bond_params_dict");
      scitbx::af::boost_python::shared_wrapper<bond_params_dict, rir>::wrap(
        "bond_params_table");
    }
  };

  struct bond_simple_proxy_wrappers
  {
    typedef bond_simple_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t, bases<bond_params> >("bond_simple_proxy", no_init)
        .def(init<af::tiny<unsigned, 2> const&, double, double>(
          (arg_("i_seqs"), arg_("distance_ideal"), arg_("weight"))))
        .def("sort_i_seqs", &w_t::sort_i_seqs)
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<bond_simple_proxy, rir>::wrap(
          "shared_bond_simple_proxy");
      }
    }
  };

  struct bond_asu_proxy_wrappers
  {
    typedef bond_asu_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<bond_params, asu_mapping_index_pair> >(
            "bond_asu_proxy", no_init)
        .def(init<asu_mapping_index_pair const&, double, double>(
          (arg_("pair"), arg_("distance_ideal"), arg_("weight"))))
        .def(init<asu_mapping_index_pair const&, bond_params const&>(
          (arg_("pair"), arg_("params"))))
        .def("as_simple_proxy", &w_t::as_simple_proxy)
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<bond_asu_proxy, rir>::wrap(
          "shared_bond_asu_proxy");
      }
    }
  };

  struct bond_wrappers
  {
    typedef bond w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t, bases<bond_params> >("bond", no_init)
        .def(init<af::tiny<scitbx::vec3<double>, 2> const&, double, double>(
          (arg_("sites"), arg_("distance_ideal"), arg_("weight"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  bond_simple_proxy const&>(
          (arg_("sites_cart"), arg_("proxy"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  asu_mappings const&,
                  bond_asu_proxy const&>(
          (arg_("sites_cart"), arg_("asu_mappings"), arg_("proxy"))))
        .add_property("sites", make_getter(&w_t::sites, rbv()))
        .def_readonly("distance_model", &w_t::distance_model)
        .def_readonly("delta", &w_t::delta)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients)
      ;
    }
  };

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    bond_residual_sum_overloads, bond_residual_sum, 3, 4)

  void
  wrap_all()
  {
    using namespace boost::python;
    bond_params_wrappers::wrap();
    bond_params_table_wrappers::wrap();
    bond_simple_proxy_wrappers::wrap();
    bond_asu_proxy_wrappers::wrap();
    bond_wrappers::wrap();
    def("extract_bond_params", extract_bond_params, (
      (arg_("n_seq"), arg_("bond_simple_proxies"))));
    def("bond_distances_model",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_simple_proxy> const&))
      bond_distances_model,
      (arg_("sites_cart"), arg_("proxies")));
    def("bond_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_simple_proxy> const&))
      bond_deltas,
      (arg_("sites_cart"), arg_("proxies")));
    def("bond_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_simple_proxy> const&))
      bond_residuals,
      (arg_("sites_cart"), arg_("proxies")));
    def("bond_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_simple_proxy> const&,
        af::ref<scitbx::vec3<double> > const&))
      bond_residual_sum,
      (arg_("sites_cart"), arg_("proxies"), arg_("gradient_array")));
    def("bond_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        bond_sorted_asu_proxies_base const&))
      bond_deltas,
      (arg_("sites_cart"), arg_("sorted_asu_proxies")));
    def("bond_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        bond_sorted_asu_proxies_base const&))
      bond_residuals,
      (arg_("sites_cart"), arg_("sorted_asu_proxies")));
    def("bond_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        bond_sorted_asu_proxies_base const&,
        af::ref<scitbx::vec3<double> > const&,
        bool)) bond_residual_sum,
      bond_residual_sum_overloads((
        arg_("sites_cart"),
        arg_("sorted_asu_proxies"),
        arg_("gradient_array"),
        arg_("disable_cache")=false)));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_bond() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
