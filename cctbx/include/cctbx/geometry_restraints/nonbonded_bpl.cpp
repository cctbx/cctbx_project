#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
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
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t>("nonbonded_params", no_init)
        .def(init<optional<double, double, double, double> >(
            (arg_("factor_1_4_interactions")=2/3.,
             arg_("const_shrink_1_4_interactions")=0,
             arg_("default_distance")=0,
             arg_("minimum_distance")=0)))
        .def_readonly("distance_table", &w_t::distance_table)
        .def_readonly("radius_table", &w_t::radius_table)
        .def_readwrite("factor_1_4_interactions",
                  &w_t::factor_1_4_interactions)
        .def_readwrite("const_shrink_1_4_interactions",
                  &w_t::const_shrink_1_4_interactions)
        .def_readwrite("default_distance", &w_t::default_distance)
        .def_readwrite("minimum_distance", &w_t::minimum_distance)
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
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("nonbonded_simple_proxy", no_init)
        .def(init<af::tiny<unsigned, 2> const&, double>(
            (arg_("i_seqs"), arg_("vdw_distance"))))
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .def_readwrite("vdw_distance", &w_t::vdw_distance)
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<w_t, rir>::wrap(
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
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t, bases<asu_mapping_index_pair> >(
            "nonbonded_asu_proxy", no_init)
        .def(init<asu_mapping_index_pair const&, double>(
          (arg_("pair"), arg_("vdw_distance"))))
        .def_readwrite("vdw_distance", &w_t::vdw_distance)
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<w_t, rir>::wrap(
          "shared_nonbonded_asu_proxy");
      }
    }
  };

  struct repulsion_function_wrappers
  {
    typedef repulsion_function w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t>("repulsion_function", no_init)
        .def(init<optional<double, double, double, double> >(
          (arg_("c_rep"), arg_("k_rep"), arg_("irexp"), arg_("rexp"))))
        .def_readonly("c_rep", &w_t::c_rep)
        .def_readonly("k_rep", &w_t::k_rep)
        .def_readonly("irexp", &w_t::irexp)
        .def_readonly("rexp", &w_t::rexp)
      ;
    }
  };

  struct nonbonded_wrappers
  {
    typedef nonbonded w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("nonbonded", no_init)
        .def(init<af::tiny<scitbx::vec3<double>, 2> const&,
                  double,
                  optional<repulsion_function const& > >(
          (arg_("sites"), arg_("vdw_distance"), arg_("function"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  nonbonded_simple_proxy const&,
                  optional<repulsion_function const& > >(
          (arg_("sites_cart"), arg_("proxy"), arg_("function"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  asu_mappings const&,
                  nonbonded_asu_proxy const&,
                  optional<repulsion_function const& > >(
          (arg_("sites_cart"), arg_("asu_mappings"), arg_("proxy"),
           arg_("function"))))
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

  struct nonbonded_sorted_asu_proxies_wrappers
  {
    typedef nonbonded_sorted_asu_proxies w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("nonbonded_sorted_asu_proxies", no_init)
        .def(init<
          boost::shared_ptr<asu_mappings> const&>(
            (arg_("asu_mappings"))))
        .def("asu_mappings", &w_t::asu_mappings, ccr())
        .def("process",
          (bool(w_t::*)(nonbonded_simple_proxy const&)) &w_t::process,
            (arg_("proxy")))
        .def("process",
          (void(w_t::*)(af::const_ref<nonbonded_simple_proxy> const&))
            &w_t::process,
          (arg_("proxies")))
        .def("process",
          (bool(w_t::*)(nonbonded_asu_proxy const&)) &w_t::process,
            (arg_("proxy")))
        .def("process",
          (void(w_t::*)(af::const_ref<nonbonded_asu_proxy> const&))
            &w_t::process,
          (arg_("proxies")))
        .def("n_total", &w_t::n_total)
        .def_readonly("simple", &w_t::simple)
        .def_readonly("asu", &w_t::asu)
      ;
    }
  };

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    nonbonded_deltas_overloads, nonbonded_deltas, 2, 3)
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    nonbonded_residuals_overloads, nonbonded_residuals, 2, 3)
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    nonbonded_residual_sum_overloads_1, nonbonded_residual_sum, 3, 4)
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    nonbonded_residual_sum_overloads_2, nonbonded_residual_sum, 3, 5)

  void
  wrap_all()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_; // gcc 2.96 workaround
    nonbonded_params_wrappers::wrap();
    nonbonded_simple_proxy_wrappers::wrap();
    nonbonded_asu_proxy_wrappers::wrap();
    repulsion_function_wrappers::wrap();
    nonbonded_wrappers::wrap();
    nonbonded_sorted_asu_proxies_wrappers::wrap();
    def("nonbonded_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<nonbonded_simple_proxy> const&,
        repulsion_function const& function)) nonbonded_deltas,
      nonbonded_deltas_overloads(
        (arg_("sites_cart"), arg_("proxies"), arg_("function"))));
    def("nonbonded_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<nonbonded_simple_proxy> const&,
        repulsion_function const& function)) nonbonded_residuals,
      nonbonded_residuals_overloads(
        (arg_("sites_cart"), arg_("proxies"), arg_("function"))));
    def("nonbonded_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<nonbonded_simple_proxy> const&,
        af::ref<scitbx::vec3<double> > const&,
        repulsion_function const& function)) nonbonded_residual_sum,
      nonbonded_residual_sum_overloads_1(
        (arg_("sites_cart"), arg_("proxies"), arg_("gradient_array"),
         arg_("function"))));
    def("nonbonded_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        nonbonded_sorted_asu_proxies const&,
        repulsion_function const& function)) nonbonded_deltas,
      nonbonded_deltas_overloads(
        (arg_("sites_cart"), arg_("sorted_asu_proxies"),
         arg_("function"))));
    def("nonbonded_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        nonbonded_sorted_asu_proxies const&,
        repulsion_function const& function)) nonbonded_residuals,
      nonbonded_residuals_overloads(
        (arg_("sites_cart"), arg_("sorted_asu_proxies"),
         arg_("function"))));
    def("nonbonded_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        nonbonded_sorted_asu_proxies const&,
        af::ref<scitbx::vec3<double> > const&,
        repulsion_function const&,
        bool)) nonbonded_residual_sum,
      nonbonded_residual_sum_overloads_2(
        (arg_("sites_cart"), arg_("sorted_asu_proxies"),
         arg_("gradient_array"), arg_("function"),
         arg_("disable_cache")=false)));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_nonbonded() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
