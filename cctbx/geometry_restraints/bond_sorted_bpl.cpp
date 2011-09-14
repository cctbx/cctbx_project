#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <cctbx/geometry_restraints/bond_sorted.h>

namespace cctbx { namespace geometry_restraints {
namespace {

  struct bond_sorted_asu_proxies_base_wrappers
  {
    typedef bond_sorted_asu_proxies_base w_t;

    static void
    process_bond_asu_proxy(
      w_t& O,
      bond_asu_proxy const& proxy)
    {
      O.process(proxy);
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("bond_sorted_asu_proxies_base", no_init)
        .def(init<
          boost::shared_ptr<asu_mappings> const&>(
            (arg("asu_mappings"))))
        .def("asu_mappings", &w_t::asu_mappings, ccr())
        .def("process", (void(w_t::*)(bond_simple_proxy const&)) &w_t::process,
          (arg("proxy")))
        .def("process",
          (void(w_t::*)(af::const_ref<bond_simple_proxy> const&))
            &w_t::process,
          (arg("proxies")))
        .def("process", process_bond_asu_proxy, (arg("proxy")))
        .def("process",
          (void(w_t::*)(af::const_ref<bond_asu_proxy> const&))
            &w_t::process,
          (arg("proxies")))
        .def("push_back",
          (void(w_t::*)(bond_asu_proxy const&)) &w_t::push_back,
            (arg("proxy")))
        .def("push_back",
          (void(w_t::*)(af::const_ref<bond_asu_proxy> const&))
            &w_t::push_back,
          (arg("proxies")))
        .def("n_total", &w_t::n_total)
        .def_readonly("simple", &w_t::simple)
        .def_readonly("asu", &w_t::asu)
      ;
    }
  };

  struct bond_sorted_asu_proxies_wrappers
  {
    typedef bond_sorted_asu_proxies w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t, bases<bond_sorted_asu_proxies_base> >(
        "bond_sorted_asu_proxies", no_init)
          .def(init<
            boost::shared_ptr<asu_mappings> const&>(
              (arg("asu_mappings"))))
          .def(init<
            af::const_ref<bond_params_dict> const&>(
              (arg("bond_params_table"))))
          .def(init<
            af::const_ref<bond_params_dict> const&,
            crystal::pair_asu_table<> const&>(
              (arg("bond_params_table"), arg("bond_asu_table"))))
          .def(init<
            crystal::pair_asu_table<> const&>(
              (arg("pair_asu_table"))))
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    bond_sorted_asu_proxies_base_wrappers::wrap();
    bond_sorted_asu_proxies_wrappers::wrap();
    def("add_pairs", add_pairs, (
      arg("pair_asu_table"), arg("bond_simple_proxies")));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_bond_sorted() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
