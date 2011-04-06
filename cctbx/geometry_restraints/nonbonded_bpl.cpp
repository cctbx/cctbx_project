#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/geometry_restraints/nonbonded.h>

namespace cctbx { namespace geometry_restraints {
namespace {

  struct nonbonded_params_wrappers
  {
    typedef nonbonded_params w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("nonbonded_params", no_init)
        .def(init<optional<double, double, double, double, double> >(
            (arg("factor_1_4_interactions")=2/3.,
             arg("const_shrink_1_4_interactions")=0,
             arg("default_distance")=0,
             arg("minimum_distance")=0,
             arg("const_shrink_donor_acceptor")=0)))
        .def("find_max_vdw_distance", &w_t::find_max_vdw_distance,
          (arg("nonbonded_types")))
        .def_readonly("distance_table", &w_t::distance_table)
        .def_readonly("radius_table", &w_t::radius_table)
        .def_readonly("donor_acceptor_table", &w_t::donor_acceptor_table)
        .def_readwrite("factor_1_4_interactions",
                  &w_t::factor_1_4_interactions)
        .def_readwrite("const_shrink_1_4_interactions",
                  &w_t::const_shrink_1_4_interactions)
        .def_readwrite("default_distance", &w_t::default_distance)
        .def_readwrite("minimum_distance", &w_t::minimum_distance)
        .def_readwrite("const_shrink_donor_acceptor",
                  &w_t::const_shrink_donor_acceptor)
      ;
    }
  };

  struct nonbonded_simple_proxy_wrappers
  {
    typedef nonbonded_simple_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("nonbonded_simple_proxy", no_init)
        .def(init<af::tiny<unsigned, 2> const&, double>(
            (arg("i_seqs"), arg("vdw_distance"))))
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .def_readwrite("vdw_distance", &w_t::vdw_distance)
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<
          nonbonded_simple_proxy, rir>::wrap(
            "shared_nonbonded_simple_proxy");
      }
    }
  };

  struct nonbonded_asu_proxy_wrappers
  {
    typedef nonbonded_asu_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<asu_mapping_index_pair> >(
            "nonbonded_asu_proxy", no_init)
        .def(init<asu_mapping_index_pair const&, double>(
          (arg("pair"), arg("vdw_distance"))))
        .def_readwrite("vdw_distance", &w_t::vdw_distance)
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<nonbonded_asu_proxy, rir>::wrap(
          "shared_nonbonded_asu_proxy");
      }
    }
  };

  struct prolsq_repulsion_function_wrappers
  {
    typedef prolsq_repulsion_function w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("prolsq_repulsion_function", no_init)
        .def(init<optional<double, double, double, double> >(
          (arg("c_rep")=16,
           arg("k_rep")=1,
           arg("irexp")=1,
           arg("rexp")=4)))
        .def_readonly("c_rep", &w_t::c_rep)
        .def_readonly("k_rep", &w_t::k_rep)
        .def_readonly("irexp", &w_t::irexp)
        .def_readonly("rexp", &w_t::rexp)
        .def("residual",
          (double (w_t::*)(double, double) const) &w_t::residual,
            (arg("vdw_distance"), arg("delta")))
      ;
    }
  };

  struct inverse_power_repulsion_function_wrappers
  {
    typedef inverse_power_repulsion_function w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("inverse_power_repulsion_function", no_init)
        .def(init<double, optional<double, double> >(
          (arg("nonbonded_distance_cutoff"),
           arg("k_rep")=1,
           arg("irexp")=1)))
        .def_readonly("nonbonded_distance_cutoff",
                 &w_t::nonbonded_distance_cutoff)
        .def_readonly("k_rep", &w_t::k_rep)
        .def_readonly("irexp", &w_t::irexp)
        .def("residual",
          (double (w_t::*)(double, double) const) &w_t::residual,
            (arg("vdw_distance"), arg("delta")))
      ;
    }
  };

  struct cos_repulsion_function_wrappers
  {
    typedef cos_repulsion_function w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("cos_repulsion_function", no_init)
        .def(init<double, optional<double> >(
          (arg("max_residual"),
           arg("exponent")=1)))
        .def_readonly("max_residual", &w_t::max_residual)
        .def_readonly("exponent", &w_t::exponent)
        .def("residual",
          (double (w_t::*)(double, double) const) &w_t::residual,
            (arg("vdw_distance"), arg("delta")))
      ;
    }
  };

  struct gaussian_repulsion_function_wrappers
  {
    typedef gaussian_repulsion_function w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("gaussian_repulsion_function", no_init)
        .def(init<double, optional<double> >(
          (arg("max_residual"),
           arg("norm_height_at_vdw_distance")=0.1)))
        .def_readonly("max_residual", &w_t::max_residual)
        .def("norm_height_at_vdw_distance", &w_t::norm_height_at_vdw_distance)
        .def("residual",
          (double (w_t::*)(double, double) const) &w_t::residual,
            (arg("vdw_distance"), arg("delta")))
      ;
    }
  };

  template <typename NonbondedFunction>
  struct nonbonded_wrappers
  {
    typedef nonbonded<NonbondedFunction> w_t;

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>(python_name, no_init)
        .def(init<af::tiny<scitbx::vec3<double>, 2> const&,
                  double,
                  NonbondedFunction const&>(
          (arg("sites"), arg("vdw_distance"), arg("function"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  nonbonded_simple_proxy const&,
                  NonbondedFunction const&>(
          (arg("sites_cart"), arg("proxy"), arg("function"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  asu_mappings const&,
                  nonbonded_asu_proxy const&,
                  NonbondedFunction const&>(
          (arg("sites_cart"), arg("asu_mappings"), arg("proxy"),
           arg("function"))))
        .add_property("sites", make_getter(&w_t::sites, rbv()))
        .def_readonly("vdw_distance", &w_t::vdw_distance)
        .def_readonly("function", &w_t::function)
        .add_property("diff_vec", make_getter(&w_t::diff_vec, rbv()))
        .def_readonly("delta", &w_t::delta)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients)
      ;
    }
  };

  template <typename NonbondedFunction>
  void
  wrap_functions(scitbx::type_holder<NonbondedFunction> const&)
  {
    using boost::python::def;
    using boost::python::arg;
    def("nonbonded_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<nonbonded_simple_proxy> const&)) nonbonded_deltas,
      (arg("sites_cart"), arg("proxies")));
    def("nonbonded_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<nonbonded_simple_proxy> const&,
        NonbondedFunction const& function)) nonbonded_residuals,
      (arg("sites_cart"), arg("proxies"), arg("function")));
    def("nonbonded_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<nonbonded_simple_proxy> const&,
        af::ref<scitbx::vec3<double> > const&,
        NonbondedFunction const& function)) nonbonded_residual_sum,
      (arg("sites_cart"), arg("proxies"), arg("gradient_array"),
       arg("function")));
    def("nonbonded_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        nonbonded_sorted_asu_proxies_base const&)) nonbonded_deltas,
      (arg("sites_cart"), arg("sorted_asu_proxies")));
    def("nonbonded_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        nonbonded_sorted_asu_proxies_base const&,
        NonbondedFunction const& function)) nonbonded_residuals,
      (arg("sites_cart"), arg("sorted_asu_proxies"), arg("function")));
    def("nonbonded_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        nonbonded_sorted_asu_proxies_base const&,
        af::ref<scitbx::vec3<double> > const&,
        NonbondedFunction const&,
        bool)) nonbonded_residual_sum, (
          arg("sites_cart"),
          arg("sorted_asu_proxies"),
          arg("gradient_array"),
          arg("function"),
          arg("disable_cache")=false));
  }

  void
  wrap_all()
  {
    nonbonded_params_wrappers::wrap();
    nonbonded_simple_proxy_wrappers::wrap();
    nonbonded_asu_proxy_wrappers::wrap();
    prolsq_repulsion_function_wrappers::wrap();
    inverse_power_repulsion_function_wrappers::wrap();
    cos_repulsion_function_wrappers::wrap();
    gaussian_repulsion_function_wrappers::wrap();
    nonbonded_wrappers<prolsq_repulsion_function>::wrap(
      "nonbonded_prolsq");
    nonbonded_wrappers<inverse_power_repulsion_function>::wrap(
      "nonbonded_inverse_power");
    nonbonded_wrappers<cos_repulsion_function>::wrap(
      "nonbonded_cos");
    nonbonded_wrappers<gaussian_repulsion_function>::wrap(
      "nonbonded_gaussian");
    wrap_functions(scitbx::type_holder<prolsq_repulsion_function>());
    wrap_functions(scitbx::type_holder<inverse_power_repulsion_function>());
    wrap_functions(scitbx::type_holder<cos_repulsion_function>());
    wrap_functions(scitbx::type_holder<gaussian_repulsion_function>());
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_nonbonded() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
