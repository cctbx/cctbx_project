#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/restraints/bond.h>

namespace cctbx { namespace restraints {
namespace {

  struct bond_proxy_wrappers
  {
    typedef bond_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("bond_proxy", no_init)
        .def(init<af::tiny<std::size_t, 2> const&, double, double>(
          (arg_("i_seqs"), arg_("distance_ideal"), arg_("weight"))))
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .def_readonly("distance_ideal", &w_t::distance_ideal)
        .def_readonly("weight", &w_t::weight)
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_bond_proxy");
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
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t>("bond_sym_proxy", no_init)
        .def(init<
          direct_space_asu::asu_mapping_index_pair const&,
          double,
          double>(
            (arg_("pair"), arg_("distance_ideal"), arg_("weight"))))
        .def_readonly("pair", &w_t::pair)
        .def_readonly("distance_ideal", &w_t::distance_ideal)
        .def_readonly("weight", &w_t::weight)
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_bond_sym_proxy");
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
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("bond", no_init)
        .def(init<af::tiny<scitbx::vec3<double>, 2> const&, double, double>(
          (arg_("sites"), arg_("distance_ideal"), arg_("weight"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  bond_proxy const&>(
          (arg_("sites_cart"), arg_("proxy"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  direct_space_asu::asu_mappings<> const&,
                  bond_sym_proxy const&>(
          (arg_("sites_cart"), arg_("asu_mappings"), arg_("proxy"))))
        .add_property("sites", make_getter(&w_t::sites, rbv()))
        .def_readonly("distance_ideal", &w_t::distance_ideal)
        .def_readonly("weight", &w_t::weight)
        .def_readonly("distance_model", &w_t::distance_model)
        .def_readonly("delta", &w_t::delta)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients)
      ;
    }
  };

  struct bond_sorted_proxies_wrappers
  {
    typedef bond_sorted_proxies w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("bond_sorted_proxies", no_init)
        .def(init<
          boost::shared_ptr<direct_space_asu::asu_mappings<> > const&>(
            (arg_("asu_mappings"))))
        .def("asu_mappings", &w_t::asu_mappings, ccr())
        .def("process", (bool(w_t::*)(bond_proxy const&)) &w_t::process,
          (arg_("proxy")))
        .def("process", (bool(w_t::*)(bond_sym_proxy const&)) &w_t::process,
          (arg_("proxy")))
        .def("n_total", &w_t::n_total)
        .def_readonly("proxies", &w_t::proxies)
        .def_readonly("sym_proxies", &w_t::sym_proxies)
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_; // gcc 2.96 workaround
    bond_proxy_wrappers::wrap();
    bond_sym_proxy_wrappers::wrap();
    bond_wrappers::wrap();
    bond_sorted_proxies_wrappers::wrap();
    def("bond_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_proxy> const&))
      bond_deltas,
      (arg_("sites_cart"), arg_("proxies")));
    def("bond_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_proxy> const&))
      bond_residuals,
      (arg_("sites_cart"), arg_("proxies")));
    def("bond_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_proxy> const&,
        af::ref<scitbx::vec3<double> > const&))
      bond_residual_sum,
      (arg_("sites_cart"), arg_("proxies"), arg_("gradient_array")));
    def("bond_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        direct_space_asu::asu_mappings<> const&,
        af::const_ref<bond_sym_proxy> const&))
      bond_deltas,
      (arg_("sites_cart"), arg_("asu_mappings"), arg_("proxies")));
    def("bond_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        direct_space_asu::asu_mappings<> const&,
        af::const_ref<bond_sym_proxy> const&))
      bond_residuals,
      (arg_("sites_cart"), arg_("asu_mappings"), arg_("proxies")));
    def("bond_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        direct_space_asu::asu_mappings<> const&,
        af::const_ref<bond_sym_proxy> const&,
        af::ref<scitbx::vec3<double> > const&))
      bond_residual_sum,
      (arg_("sites_cart"),
       arg_("asu_mappings"),
       arg_("proxies"),
       arg_("gradient_array")));
    def("bond_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        bond_sorted_proxies const&,
        af::ref<scitbx::vec3<double> > const&))
      bond_residual_sum,
      (arg_("sites_cart"),
       arg_("sorted_proxies"),
       arg_("gradient_array")));
    def("bond_sets", bond_sets,
      (arg_("n_sites"), arg_("bond_proxies")));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_bond() { wrap_all(); }

}}} // namespace cctbx::restraints::boost_python
