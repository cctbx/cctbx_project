#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <cctbx/geometry_restraints/bond_sorted.h>
#include <cctbx/geometry_restraints/sorted_asu_proxies.h>

namespace cctbx { namespace geometry_restraints {
namespace {

  struct bond_sorted_asu_proxies_base_wrappers : boost::python::pickle_suite
  {
    typedef bond_sorted_asu_proxies_base w_t;

    static boost::python::tuple
      getstate(w_t const& self)
    {
      return boost::python::make_tuple(
        self.asu_mappings(),
        self.simple,
        self.asu
        );
    }

    static void
      setstate(w_t& self, boost::python::tuple state)
    {
      self.asu_mappings_owner_ = boost::python::extract< boost::shared_ptr<direct_space_asu::asu_mappings<> > >(state[0]);
      self.simple = boost::python::extract< af::shared<bond_simple_proxy> >(state[1]);
      self.asu = boost::python::extract< af::shared<bond_asu_proxy> >(state[2]);
    }

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
        .def_pickle(bond_sorted_asu_proxies_base_wrappers())
        ;
    }
  };

  struct bond_sorted_asu_proxies_wrappers : boost::python::pickle_suite
  {
    typedef bond_sorted_asu_proxies w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
        return boost::python::make_tuple(self.asu_mappings());
    }

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
          .def_pickle(bond_sorted_asu_proxies_wrappers())
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
