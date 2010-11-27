#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/array_family/selections.h>
#include <scitbx/stl/map_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <cctbx/geometry_restraints/bond_misc.h>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(
  cctbx::geometry_restraints::bond_asu_proxy)

namespace scitbx { namespace af { namespace boost_python {

  using cctbx::geometry_restraints::bond_simple_proxy;

  template <>
  struct shared_wrapper_default_element<bond_simple_proxy>
  {
    static bond_simple_proxy
    get()
    {
      return bond_simple_proxy(
        af::tiny<unsigned, 2>(0, 0), 0., 0., 0.);
    }
  };

}}}

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
        .def(init<double, double, optional<double> >((
          arg("distance_ideal"), arg("weight"), arg("slack")=0)))
        .def("scale_weight", &w_t::scale_weight, (arg("factor")))
        .def_readwrite("distance_ideal", &w_t::distance_ideal)
        .def_readwrite("weight", &w_t::weight)
        .def_readwrite("slack", &w_t::slack)
      ;
    }
  };

  struct bond_params_table_wrappers
  {
    static void
    update(
      bond_params_table& self,
      unsigned i_seq,
      unsigned j_seq,
      bond_params const& params)
    {
      if (i_seq <= j_seq) self[i_seq][j_seq] = params;
      else                self[j_seq][i_seq] = params;
    }

    static double
    mean_residual(
      af::const_ref<bond_params_dict> const& self,
      double bond_stretch_factor)
    {
      double sum = 0;
      unsigned n = 0;
      for(unsigned i_seq=0;i_seq<self.size();i_seq++) {
        for(bond_params_dict::const_iterator
              dict_i=self[i_seq].begin();
              dict_i!=self[i_seq].end();
              dict_i++) {
          bond_params const& params = dict_i->second;
          double delta = params.distance_ideal * bond_stretch_factor;
          sum += params.weight * delta * delta;
        }
        n += self[i_seq].size();
      }
      if (n == 0) return 0;
      return sum / static_cast<double>(n);
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      scitbx::stl::boost_python::map_wrapper<bond_params_dict, rir>::wrap(
        "bond_params_dict");
      scitbx::af::boost_python::shared_wrapper<bond_params_dict, rir>::wrap(
        "bond_params_table")
        .def("update", update, (
          arg("self"),
          arg("i_seq"),
          arg("j_seq"),
          arg("params")))
        .def("mean_residual", mean_residual, (
          arg("self"),
          arg("bond_stretch_factor")))
        .def("proxy_select",
          (bond_params_table(*)(
            af::const_ref<bond_params_dict> const&,
            af::const_ref<std::size_t> const&))
              scitbx::af::array_of_map_proxy_select, (
          arg("iselection")))
        .def("proxy_remove",
          (bond_params_table(*)(
            af::const_ref<bond_params_dict> const&,
            af::const_ref<bool> const&))
              scitbx::af::array_of_map_proxy_remove, (
          arg("selection")))
      ;
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
        .def(init<
          af::tiny<unsigned, 2> const&, double, double, optional<double> >((
            arg("i_seqs"),
            arg("distance_ideal"),
            arg("weight"),
            arg("slack")=0)))
        .def(init<
          af::tiny<unsigned, 2> const&,
          sgtbx::rt_mx const&,
          double,
          double,
          optional<double> >((
            arg("i_seqs"),
            arg("rt_mx_ji"),
            arg("distance_ideal"),
            arg("weight"),
            arg("slack")=0)))
        .def("sort_i_seqs", &w_t::sort_i_seqs)
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .add_property("rt_mx_ji", make_getter(&w_t::rt_mx_ji, rbv()))
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<bond_simple_proxy, rir>::wrap(
          "shared_bond_simple_proxy");
      }
    }
  };

  struct bond_sym_proxy_wrappers
  {
    typedef bond_sym_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t, bases<bond_params> >("bond_sym_proxy", no_init)
        .def(init<
          af::tiny<unsigned, 2> const&,
          sgtbx::rt_mx const&,
          double,
          double,
          optional<double> >((
            arg("i_seqs"),
            arg("rt_mx_ji"),
            arg("distance_ideal"),
            arg("weight"),
            arg("slack")=0)))
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .def_readonly("rt_mx_ji", &w_t::rt_mx_ji)
      ;
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
        .def(init<
          asu_mapping_index_pair const&, double, double, optional<double> >((
            arg("pair"),
            arg("distance_ideal"),
            arg("weight"),
            arg("slack")=0)))
        .def(init<asu_mapping_index_pair const&, bond_params const&>(
          (arg("pair"), arg("params"))))
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
      class_<w_t, bases<bond_params, cctbx::geometry::distance<double> > >(
        "bond", no_init)
        .def(init<
          af::tiny<scitbx::vec3<double>, 2> const&,
          double,
          double,
          optional<double> >((
            arg("sites"),
            arg("distance_ideal"),
            arg("weight"),
            arg("slack")=0)))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  bond_simple_proxy const&>(
          (arg("sites_cart"), arg("proxy"))))
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<scitbx::vec3<double> > const&,
                  bond_simple_proxy const&>(
          (arg("unit_cell"), arg("sites_cart"), arg("proxy"))))
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<scitbx::vec3<double> > const&,
                  bond_sym_proxy const&>(
          (arg("unit_cell"), arg("sites_cart"), arg("proxy"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  asu_mappings const&,
                  bond_asu_proxy const&>(
          (arg("sites_cart"), arg("asu_mappings"), arg("proxy"))))
        .add_property("sites", make_getter(&w_t::sites, rbv()))
        .def_readonly("distance_model", &w_t::distance_model)
        .def_readonly("delta", &w_t::delta)
        .def_readonly("delta_slack", &w_t::delta_slack)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients)
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    bond_params_wrappers::wrap();
    bond_params_table_wrappers::wrap();
    bond_simple_proxy_wrappers::wrap();
    bond_sym_proxy_wrappers::wrap();
    bond_asu_proxy_wrappers::wrap();
    bond_wrappers::wrap();
    def("extract_bond_params", extract_bond_params, (
      (arg("n_seq"), arg("bond_simple_proxies"))));
    def("bond_distances_model",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_simple_proxy> const&))
      bond_distances_model,
      (arg("sites_cart"), arg("proxies")));
    def("bond_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_simple_proxy> const&))
      bond_deltas,
      (arg("sites_cart"), arg("proxies")));
    def("bond_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_simple_proxy> const&))
      bond_residuals,
      (arg("sites_cart"), arg("proxies")));
    def("bond_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_simple_proxy> const&,
        af::ref<scitbx::vec3<double> > const&))
      bond_residual_sum,
      (arg("sites_cart"), arg("proxies"), arg("gradient_array")));
    def("bond_deltas",
      (af::shared<double>(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_simple_proxy> const&))
      bond_deltas,
      (arg("unit_cell"), arg("sites_cart"), arg("proxies")));
    def("bond_residuals",
      (af::shared<double>(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_simple_proxy> const&))
      bond_residuals,
      (arg("unit_cell"), arg("sites_cart"), arg("proxies")));
    def("bond_distances_model",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        bond_sorted_asu_proxies_base const&))
      bond_distances_model,
      (arg("sites_cart"), arg("sorted_asu_proxies")));
    def("bond_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        bond_sorted_asu_proxies_base const&))
      bond_deltas,
      (arg("sites_cart"), arg("sorted_asu_proxies")));
    def("bond_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        bond_sorted_asu_proxies_base const&))
      bond_residuals,
      (arg("sites_cart"), arg("sorted_asu_proxies")));
    def("bond_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        bond_sorted_asu_proxies_base const&,
        af::ref<scitbx::vec3<double> > const&,
        bool)) bond_residual_sum, (
          arg("sites_cart"),
          arg("sorted_asu_proxies"),
          arg("gradient_array"),
          arg("disable_cache")=false));

    def("home_restraints_summation_skip_special_positions",
      home_restraints_summation_skip_special_positions, (
        arg("sites_cart"),
        arg("gradients"),
        arg("site_symmetry_table_indices"),
        arg("home_sites_cart"),
        arg("iselection"),
        arg("weight"),
        arg("slack")));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_bond() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
